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
        <table style="width:100%">
    """)
    footer = textwrap.dedent("""
        </table>
      </body>
    </html>
    """)
    lines = ["      " + line for line in lines]
    text = header + "\n".join(lines) + footer
    Path('index.html').write_text(text)


def link(url, text):
    return '<a href="' + url + '">' + text + '</a>'


def td(txt):
    return "<td>" + txt + "</td>"


def th(txt):
    return "<th>" + txt + "</th>"


def tr(entries):
    return "<tr>" + "  ".join(entries) + "</tr>"


def generate_summary(output_dir):

    lines = []

    table_header = ["status", "recipe"]
    lines.append(tr(th(txt) for txt in table_header))

    for recipe_dir in sorted(Path(output_dir).glob('recipe_*')):
        log = recipe_dir / 'run' / 'main_log.txt'
        success = log.read_text().endswith('Run was successful\n')

        entry = []
        entry.append('success' if success else 'failed')
        entry.append(link(recipe_dir.name, recipe_dir.name))

        entry_txt = tr(td(txt) for txt in entry)
        lines.append(entry_txt)

    return lines


def main(output_dir):

    write_summary(generate_summary(output_dir))


if __name__ == '__main__':
    main(sys.argv[1])
