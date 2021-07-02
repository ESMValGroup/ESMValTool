"""Write an index.html file in a directory containing recipe runs."""
import argparse
import textwrap
from pathlib import Path


def write_index_html(lines, output_dir):
    """Write lines to index.html."""
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

    index_file = output_dir / 'index.html'
    index_file.write_text(text)
    print(f"Wrote file://{index_file.absolute()}")


def link(url, text):
    """Format text as html link."""
    return '<a href="' + url + '">' + text + '</a>'


def tr(entries):
    """Format text entries as html table row."""
    return "<tr>" + "  ".join(entries) + "</tr>"


def th(txt):
    """Format text as html table header."""
    return "<th>" + txt + "</th>"


def td(txt):
    """Format text as html table data."""
    return "<td>" + txt + "</td>"


def generate_summary(output_dir):
    """Generate the lines of text for the summary view."""
    lines = []

    column_titles = ["status", "recipe"]
    lines.append(tr(th(txt) for txt in column_titles))

    for recipe_dir in sorted(Path(output_dir).glob('recipe_*')):
        log = recipe_dir / 'run' / 'main_log.txt'
        success = log.read_text().endswith('Run was successful\n')

        entry = []
        entry.append('success' if success else 'failed')
        entry.append(link(recipe_dir.name, recipe_dir.name))

        entry_txt = tr(td(txt) for txt in entry)
        lines.append(entry_txt)

    return lines


def main():
    """Run the program."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir',
                        default='.',
                        type=Path,
                        help='ESMValTool output directory.')
    args = parser.parse_args()

    write_index_html(generate_summary(args.output_dir), args.output_dir)


if __name__ == '__main__':
    main()
