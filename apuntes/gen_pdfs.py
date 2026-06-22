#!/usr/bin/env python3
"""Pre-process markdown files and generate PDFs via pandoc+typst."""
import re, subprocess, tempfile, sys, os

def preprocess(text):
    # 1. Remove TOC section (## Índice ... ---)
    text = re.sub(
        r'^## Índice\s*\n(.*?\n)*?---',
        '---',
        text,
        count=1,
        flags=re.MULTILINE,
    )

    # 2. Fix \empty -> \emptyset
    text = text.replace(r'\empty', r'\emptyset')

    # 3. Fix multi-line inline math with \\ (convert to display aligned)
    # Match $...\\...$ spanning multiple lines
    def fix_multiline_inline_math(m):
        content = m.group(1)
        if '\\\\' in content:
            # Convert to display math with aligned
            content = content.strip()
            # Wrap in aligned environment
            return f'$$\\begin{{aligned}}\n{content}\n\\end{{aligned}}$$'
        return m.group(0)

    # Match inline math that spans multiple lines (contains actual newlines)
    text = re.sub(
        r'(?<!\$)\$(?!\$)((?:[^\$]|\n)+?\\\\(?:[^\$]|\n)*?)\$(?!\$)',
        fix_multiline_inline_math,
        text,
    )

    # 4. Fix \argmin (not standard in pandoc's tex parser)
    text = text.replace(r'\argmin', r'\operatorname{argmin}')

    # 5. Fix \iff at start of display math block
    text = re.sub(
        r'^\$\\iff',
        r'$\\Leftrightarrow',
        text,
        flags=re.MULTILINE,
    )

    # 6. Fix display math blocks with \\ that need aligned environment
    def fix_display_math(m):
        content = m.group(1)
        # If it contains \\ but no \begin{...} environment, wrap in aligned
        if '\\\\' in content and '\\begin{' not in content:
            content = content.strip()
            return f'$$\\begin{{aligned}}\n{content}\n\\end{{aligned}}$$'
        return m.group(0)

    text = re.sub(r'\$\$(.*?)\$\$', fix_display_math, text, flags=re.DOTALL)

    # 7. Replace \Vmatrix with \lVert ... \rVert pmatrix (typst compat)
    text = text.replace(r'\begin{Vmatrix}', r'\left\lVert \begin{pmatrix}')
    text = text.replace(r'\end{Vmatrix}', r'\end{pmatrix} \right\rVert')

    return text


FILES = [
    'resumen-1er-parcial.md',
    'resumen-2do-parcial.md',
    'resumen-2do-parcial-gemini.md',
]

base_dir = os.path.dirname(os.path.abspath(__file__))

for fname in FILES:
    src = os.path.join(base_dir, fname)
    pdf_name = fname.replace('.md', '.pdf')
    pdf_path = os.path.join(base_dir, pdf_name)

    print(f'Processing {fname}...')
    with open(src, 'r') as f:
        text = f.read()

    processed = preprocess(text)

    # Write to temp file
    with tempfile.NamedTemporaryFile(
        mode='w', suffix='.md', delete=False, dir=base_dir
    ) as tmp:
        tmp.write(processed)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [
                'pandoc', tmp_path,
                '-o', pdf_path,
                '--pdf-engine=typst',
                '-V', 'mainfont=New Computer Modern',
                '--toc',
                '-V', 'toc-title=Índice',
            ],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            print(f'  WARN (exit {result.returncode}):')
            # Print only errors, not warnings
            for line in result.stderr.splitlines():
                if 'error' in line.lower():
                    print(f'    {line}')
            if result.returncode >= 43:
                print(f'  FULL STDERR for debugging:')
                print(result.stderr[-2000:])
        else:
            print(f'  -> {pdf_name} OK')
            if result.stderr:
                warns = [l for l in result.stderr.splitlines() if 'WARNING' in l]
                if warns:
                    print(f'  ({len(warns)} warnings)')
    finally:
        os.unlink(tmp_path)
