import sys
import textwrap
from pathlib import Path


def write_summary(lines):

    header = textwrap.dedent("""
    <!doctype html>
    <html>
      <head>
        <title>ESMValTool recipes</title>
      </head>
      <body>
    """)
    footer = textwrap.dedent("""
      </body>
    </html>
    """)
    lines = ["    " + line for line in lines]
    text = header + "\n".join(lines) + footer
    Path('index.html').write_text(text)


def link(url, text):
    return '<a href="' + url + '">' + text + '</a>'


def generate_summary(output_dir):

    lines = []

    for output_dir in sorted(Path(output_dir).glob('recipe_*')):
        log = output_dir / 'run' / 'main_log.txt'
        success = log.read_text().endswith('Run was successful\n')

        entry = []
        entry.append('success' if success else 'failed')
        entry.append(link(output_dir.name, output_dir.name))
        entry_txt = "\t".join(entry)

        lines.append(entry_txt)

    return lines


def main(output_dir):

    write_summary(generate_summary(output_dir))


if __name__ == '__main__':
    main(sys.argv[1])
