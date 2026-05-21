from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import zipfile
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from lxml import etree


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
M = "http://schemas.openxmlformats.org/officeDocument/2006/math"
NS = {"w": W, "m": M}

INPUT = Path("_MConverter.eu_chapter02.docx")
OUTPUT = Path("chapter02_formula_fixed.docx")
MARKDOWN = Path("chapter02.md")


@dataclass(frozen=True)
class FormulaMeta:
    kind: str
    has_left_text: bool
    has_right_text: bool


@dataclass(frozen=True)
class FormulaTemplate:
    payload: etree._Element
    display: bool


@dataclass(frozen=True)
class RawFormulaGroup:
    start_paragraph: etree._Element
    end_paragraph: etree._Element
    start_text: etree._Element
    end_text: etree._Element
    start_offset: int
    end_offset: int
    preceding_math: int
    latex: str


def qn(namespace: str, tag: str) -> str:
    return f"{{{namespace}}}{tag}"


def word_text(paragraph: etree._Element) -> str:
    return "".join(paragraph.xpath(".//w:t/text()", namespaces=NS)).strip()


def has_math(paragraph: etree._Element) -> bool:
    return bool(paragraph.xpath(".//m:oMath | .//m:oMathPara", namespaces=NS))


def math_only(paragraph: etree._Element) -> bool:
    return has_math(paragraph) and not word_text(paragraph)


def paragraph_style(paragraph: etree._Element) -> str | None:
    styles = paragraph.xpath("./w:pPr/w:pStyle/@w:val", namespaces=NS)
    return styles[0] if styles else None


def markdown_formula_metadata(markdown: Path) -> list[FormulaMeta]:
    text = markdown.read_text(encoding="utf-8")
    metadata: list[FormulaMeta] = []
    for match in re.finditer(r"\$\$(.*?)\$\$", text, flags=re.S):
        line_start = text.rfind("\n", 0, match.start()) + 1
        line_end = text.find("\n", match.end())
        line_end = len(text) if line_end < 0 else line_end
        left = text[line_start : match.start()].strip()
        right = text[match.end() : line_end].strip()
        kind = "display" if not left and not right else "inline"
        metadata.append(FormulaMeta(kind, bool(left), bool(right)))
    return metadata


def pandoc_formula_templates(markdown: Path, base: Path) -> list[FormulaTemplate]:
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        raise RuntimeError("pandoc is required to recover MConverter raw TeX formulas")

    with tempfile.TemporaryDirectory(prefix="chapter02-pandoc-", dir=base) as temp_dir:
        reference = Path(temp_dir) / "formula_reference.docx"
        subprocess.run(
            [
                pandoc,
                str(markdown),
                "--from=markdown+tex_math_dollars",
                "--to=docx",
                f"--resource-path={markdown.parent}",
                f"--output={reference}",
            ],
            cwd=base,
            check=True,
            capture_output=True,
            text=True,
        )
        with zipfile.ZipFile(reference) as zref:
            reference_document = etree.fromstring(zref.read("word/document.xml"))

    templates: list[FormulaTemplate] = []
    for math in reference_document.xpath("//m:oMath", namespaces=NS):
        parent = math.getparent()
        display = parent is not None and parent.tag == qn(M, "oMathPara")
        payload = parent if display else math
        templates.append(FormulaTemplate(deepcopy(payload), display))
    return templates


def assign_formula_metadata(
    document: etree._Element,
    metadata: list[FormulaMeta],
    templates: list[FormulaTemplate],
) -> dict[etree._Element, FormulaMeta]:
    by_paragraph: dict[etree._Element, FormulaMeta] = {}
    cursor = 0
    for paragraph in document.xpath("//w:p", namespaces=NS):
        math_nodes = paragraph.xpath(".//m:oMath", namespaces=NS)
        if not math_nodes:
            continue
        if cursor < len(metadata):
            by_paragraph[paragraph] = metadata[cursor]
        cursor += len(math_nodes)
    print(f"markdown_formula_metadata={len(metadata)}")
    print(f"pandoc_formula_templates={len(templates)}")
    print(f"document_math_nodes={cursor}")
    return by_paragraph


def raw_formula_groups(document: etree._Element) -> list[RawFormulaGroup]:
    groups: list[RawFormulaGroup] = []
    opening: tuple[etree._Element, etree._Element, int, int] | None = None
    latex: list[str] = []
    needle = "$$"

    for paragraph in document.xpath("//w:p", namespaces=NS):
        for text_node in paragraph.xpath(".//w:t", namespaces=NS):
            text = text_node.text or ""
            cursor = 0
            while True:
                delimiter = text.find(needle, cursor)
                if delimiter < 0:
                    if opening is not None:
                        latex.append(text[cursor:])
                    break

                if opening is None:
                    opening = (
                        paragraph,
                        text_node,
                        delimiter,
                        len(text_node.xpath("preceding::m:oMath", namespaces=NS)),
                    )
                    cursor = delimiter + len(needle)
                    continue

                latex.append(text[cursor:delimiter])
                start_paragraph, start_text, start_offset, preceding_math = opening
                groups.append(
                    RawFormulaGroup(
                        start_paragraph,
                        paragraph,
                        start_text,
                        text_node,
                        start_offset,
                        delimiter,
                        preceding_math,
                        "".join(latex),
                    )
                )
                opening = None
                latex = []
                cursor = delimiter + len(needle)

    if opening is not None:
        raise RuntimeError("unclosed raw TeX formula remains in document.xml")
    return groups


def clear_standalone_math_delimiters(document: etree._Element) -> int:
    cleared = 0
    for paragraph in document.xpath("//w:p", namespaces=NS):
        text = "".join(paragraph.xpath(".//w:t/text()", namespaces=NS)).strip()
        if text != "$$" or has_math(paragraph):
            continue
        for text_node in paragraph.xpath(".//w:t", namespaces=NS):
            text_node.text = ""
        cleared += 1
    return cleared


def paragraph_has_visible_content(paragraph: etree._Element) -> bool:
    return bool(
        word_text(paragraph)
        or paragraph.xpath(".//m:oMath | .//m:oMathPara", namespaces=NS)
        or paragraph.xpath(".//w:drawing | .//w:pict", namespaces=NS)
    )


def clear_from_opening_delimiter(group: RawFormulaGroup) -> None:
    clear = False
    for text_node in group.start_paragraph.xpath(".//w:t", namespaces=NS):
        if text_node is group.start_text:
            text_node.text = (text_node.text or "")[: group.start_offset]
            clear = True
        elif clear:
            text_node.text = ""


def clear_through_closing_delimiter(group: RawFormulaGroup) -> None:
    clear = True
    for text_node in group.end_paragraph.xpath(".//w:t", namespaces=NS):
        if text_node is group.end_text:
            text_node.text = (text_node.text or "")[group.end_offset + 2 :]
            clear = False
        elif clear:
            text_node.text = ""


def new_formula_paragraph(template: FormulaTemplate) -> etree._Element:
    paragraph = etree.Element(qn(W, "p"))
    if template.display:
        paragraph.append(deepcopy(template.payload))
        return paragraph

    math_para = etree.Element(qn(M, "oMathPara"))
    math_para.append(deepcopy(template.payload))
    paragraph.append(math_para)
    return paragraph


def run_ancestor(element: etree._Element) -> etree._Element:
    ancestors = element.xpath("ancestor::w:r[1]", namespaces=NS)
    if not ancestors:
        raise RuntimeError("raw TeX delimiter is not inside a Word run")
    return ancestors[0]


def inline_math_payload(template: FormulaTemplate) -> etree._Element:
    if not template.display:
        return deepcopy(template.payload)
    math = template.payload.find("m:oMath", namespaces=NS)
    if math is None:
        raise RuntimeError("display formula template does not contain m:oMath")
    return deepcopy(math)


def suffix_run(run: etree._Element, text: str) -> etree._Element:
    new_run = etree.Element(qn(W, "r"))
    run_properties = run.find("w:rPr", namespaces=NS)
    if run_properties is not None:
        new_run.append(deepcopy(run_properties))
    word_text_node = etree.SubElement(new_run, qn(W, "t"))
    if text[:1].isspace() or text[-1:].isspace():
        word_text_node.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    word_text_node.text = text
    return new_run


def replace_inline_raw_group(group: RawFormulaGroup, template: FormulaTemplate) -> None:
    paragraph = group.start_paragraph
    if group.end_paragraph is not paragraph:
        raise RuntimeError("inline raw TeX formula unexpectedly spans paragraphs")

    start_run = run_ancestor(group.start_text)
    parent = start_run.getparent()
    if parent is None:
        raise RuntimeError("Word run parent is missing")

    if group.start_text is group.end_text:
        original = group.start_text.text or ""
        prefix = original[: group.start_offset]
        suffix = original[group.end_offset + 2 :]
        group.start_text.text = prefix
        insertion = parent.index(start_run) + 1
        parent.insert(insertion, inline_math_payload(template))
        if suffix:
            parent.insert(insertion + 1, suffix_run(start_run, suffix))
        return

    clear_from_opening_delimiter(group)
    clear_through_closing_delimiter(group)
    text_nodes = paragraph.xpath(".//w:t", namespaces=NS)
    clear = False
    for text_node in text_nodes:
        if text_node is group.start_text:
            clear = True
            continue
        if text_node is group.end_text:
            break
        if clear:
            text_node.text = ""
    parent.insert(parent.index(start_run) + 1, inline_math_payload(template))


def replace_display_raw_group(group: RawFormulaGroup, template: FormulaTemplate) -> None:
    start_parent = group.start_paragraph.getparent()
    end_parent = group.end_paragraph.getparent()
    if start_parent is None or start_parent is not end_parent:
        raise RuntimeError("display raw TeX formula crosses unsupported containers")

    insert_at = start_parent.index(group.start_paragraph)
    clear_from_opening_delimiter(group)
    clear_through_closing_delimiter(group)

    between = list(start_parent)[
        start_parent.index(group.start_paragraph) + 1 : start_parent.index(group.end_paragraph)
    ]
    for paragraph in between:
        start_parent.remove(paragraph)

    if paragraph_has_visible_content(group.start_paragraph):
        insert_at = start_parent.index(group.start_paragraph) + 1
    else:
        start_parent.remove(group.start_paragraph)

    if group.end_paragraph is not group.start_paragraph:
        if not paragraph_has_visible_content(group.end_paragraph):
            start_parent.remove(group.end_paragraph)

    start_parent.insert(insert_at, new_formula_paragraph(template))


def recover_raw_latex_formulas(
    document: etree._Element,
    templates: list[FormulaTemplate],
) -> int:
    groups = raw_formula_groups(document)
    for inserted, group in enumerate(groups):
        formula_index = group.preceding_math + inserted
        if formula_index >= len(templates):
            preview = " ".join(group.latex.split())[:120]
            raise RuntimeError(f"missing pandoc formula template for raw TeX: {preview}")
        template = templates[formula_index]
        if template.display or group.start_paragraph is not group.end_paragraph:
            replace_display_raw_group(group, template)
        else:
            replace_inline_raw_group(group, template)
    return len(groups)


def ensure_page_geometry(document: etree._Element) -> None:
    sect_prs = document.xpath("//w:sectPr", namespaces=NS)
    if not sect_prs:
        body = document.find("w:body", namespaces=NS)
        sect_prs = [etree.SubElement(body, qn(W, "sectPr"))]

    for sect_pr in sect_prs:
        if sect_pr.find("w:pgSz", namespaces=NS) is None:
            pg_sz = etree.Element(qn(W, "pgSz"))
            pg_sz.set(qn(W, "w"), "12240")
            pg_sz.set(qn(W, "h"), "15840")
            sect_pr.insert(0, pg_sz)
        if sect_pr.find("w:pgMar", namespaces=NS) is None:
            pg_mar = etree.Element(qn(W, "pgMar"))
            for name, value in {
                "top": "1440",
                "right": "1440",
                "bottom": "1440",
                "left": "1440",
                "header": "720",
                "footer": "720",
                "gutter": "0",
            }.items():
                pg_mar.set(qn(W, name), value)
            sect_pr.insert(1, pg_mar)


def direct_children_after_properties(paragraph: etree._Element) -> list[etree._Element]:
    return [
        deepcopy(child)
        for child in paragraph
        if etree.QName(child).localname != "pPr"
    ]


def append_paragraph_content(target: etree._Element, source: etree._Element) -> None:
    for child in direct_children_after_properties(source):
        target.append(child)


def inline_math_children(paragraph: etree._Element) -> list[etree._Element]:
    return [
        deepcopy(math)
        for math_para in paragraph.xpath("./m:oMathPara", namespaces=NS)
        for math in math_para.xpath("./m:oMath", namespaces=NS)
    ]


def replace_display_wrapper_with_inline_math(paragraph: etree._Element) -> None:
    math_children = inline_math_children(paragraph)
    for math_para in paragraph.xpath("./m:oMathPara", namespaces=NS):
        paragraph.remove(math_para)
    for math in math_children:
        paragraph.append(math)


def merge_next_text_fragment(
    container: etree._Element,
    siblings: list[etree._Element],
    index: int,
    target: etree._Element,
    metadata: dict[etree._Element, FormulaMeta],
) -> bool:
    if index + 1 >= len(siblings):
        return False
    nxt = siblings[index + 1]
    if nxt in metadata and metadata[nxt].kind == "display":
        return False
    if paragraph_style(nxt) and paragraph_style(nxt).startswith("Heading"):
        return False
    if math_only(nxt):
        return False

    append_paragraph_content(target, nxt)
    container.remove(nxt)
    siblings.pop(index + 1)
    return True


def repair_container(
    container: etree._Element,
    metadata: dict[etree._Element, FormulaMeta],
) -> tuple[int, int, int]:
    merged_paragraphs = 0
    merged_text_fragments = 0
    rewritten_paragraphs = 0
    paragraphs = list(container.findall("w:p", namespaces=NS))
    i = 0
    while i < len(paragraphs):
        current = paragraphs[i]
        meta = metadata.get(current)
        if (
            meta is not None
            and meta.kind == "inline"
            and not math_only(current)
            and current.xpath("./m:oMathPara", namespaces=NS)
        ):
            replace_display_wrapper_with_inline_math(current)
            rewritten_paragraphs += 1
            i += 1
            continue
        if not math_only(current) or meta is None or meta.kind != "inline":
            i += 1
            continue

        prev = paragraphs[i - 1] if i else None
        if meta.has_left_text and prev is not None and not (
            paragraph_style(prev) and paragraph_style(prev).startswith("Heading")
        ):
            for math in inline_math_children(current):
                prev.append(math)
            container.remove(current)
            paragraphs.pop(i)
            merged_paragraphs += 1
            previous_index = i - 1
            if meta.has_right_text and merge_next_text_fragment(
                container, paragraphs, previous_index, prev, metadata
            ):
                merged_text_fragments += 1
            continue

        replace_display_wrapper_with_inline_math(current)
        rewritten_paragraphs += 1
        if meta.has_right_text and merge_next_text_fragment(
            container, paragraphs, i, current, metadata
        ):
            merged_text_fragments += 1
        i += 1

    return merged_paragraphs, merged_text_fragments, rewritten_paragraphs


def merge_split_inline_math(
    document: etree._Element,
    metadata: dict[etree._Element, FormulaMeta],
) -> tuple[int, int, int]:
    totals = [0, 0, 0]
    containers = [document.find("w:body", namespaces=NS)]
    containers.extend(document.xpath("//w:tc", namespaces=NS))
    for container in containers:
        if container is None:
            continue
        results = repair_container(container, metadata)
        totals = [total + result for total, result in zip(totals, results)]
    return tuple(totals)


def patch_docx(source: Path, markdown: Path, output: Path) -> None:
    if output.exists():
        output.unlink()

    with zipfile.ZipFile(source) as zin:
        document = etree.fromstring(zin.read("word/document.xml"))
        templates = pandoc_formula_templates(markdown, output.parent)
        cleared_delimiters = clear_standalone_math_delimiters(document)
        recovered_raw = recover_raw_latex_formulas(document, templates)
        formula_metadata = assign_formula_metadata(
            document,
            markdown_formula_metadata(markdown),
            templates,
        )
        ensure_page_geometry(document)
        merged, text_fragments, rewritten = merge_split_inline_math(
            document,
            formula_metadata,
        )
        patched_document = etree.tostring(
            document,
            xml_declaration=True,
            encoding="UTF-8",
            standalone="yes",
        )

        with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                payload = patched_document if item.filename == "word/document.xml" else zin.read(item.filename)
                zout.writestr(item, payload)

    print(f"input={source}")
    print(f"output={output}")
    print(f"cleared_standalone_math_delimiters={cleared_delimiters}")
    print(f"recovered_raw_latex_formulas={recovered_raw}")
    print(f"merged_math_paragraphs={merged}")
    print(f"merged_text_fragments={text_fragments}")
    print(f"paragraph_start_inline_math={rewritten}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    patch_docx(base / INPUT, base / MARKDOWN, base / OUTPUT)
