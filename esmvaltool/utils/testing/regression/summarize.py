"""Write an index.html file in a directory containing recipe runs."""
import argparse
import datetime
import textwrap
from pathlib import Path


def read_resource_usage_file(recipe_dir):
    """Read resource usage from the log."""
    resource_file = recipe_dir / 'run' / 'resource_usage.txt'
    usage = {}

    if not resource_file.exists():
        return usage

    text = resource_file.read_text().strip()
    if not text:
        return usage

    lines = text.split('\n')
    for name in lines[0].split('\t'):
        usage[name] = []

    for line in lines[1:]:
        for key, value in zip(usage, line.split('\t')):
            if key != 'Date and time (UTC)':
                value = float(value)
            usage[key].append(value)

    return usage


def get_runtime_from_debug(recipe_dir):
    """Try to read the runtime from the debug log."""
    debug_file = recipe_dir / 'run' / 'main_log_debug.txt'
    if not debug_file.exists():
        return None

    text = debug_file.read_text().strip()
    if not text:
        return None

    lines = text.split('\n')
    fmt = "%Y-%m-%d %H:%M:%S"
    end_date = None
    for line in lines[::-1]:
        try:
            end_date = datetime.datetime.strptime(line[:19], fmt)
        except ValueError:
            pass
        else:
            break
    if end_date is None:
        return None

    start_date = datetime.datetime.strptime(lines[0][:19], fmt)
    runtime = end_date - start_date
    runtime = datetime.timedelta(seconds=round(runtime.total_seconds()))
    return runtime


def get_resource_usage(recipe_dir):
    """Get recipe runtime (minutes), max memory (GB), avg CPU."""
    resource_usage = read_resource_usage_file(recipe_dir)

    if not resource_usage or not resource_usage['Real time (s)']:
        runtime = get_runtime_from_debug(recipe_dir)
        runtime = "" if runtime is None else f"{runtime}"
        return [runtime, '', '']

    runtime = resource_usage['Real time (s)'][-1]
    avg_cpu = resource_usage['CPU time (s)'][-1] / runtime * 100.
    runtime = datetime.timedelta(seconds=round(runtime))
    memory = max(resource_usage['Memory (GB)'])

    return [f"{runtime}", f"{memory:.1f}", f"{avg_cpu:.1f}"]


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

    column_titles = [
        "status",
        "recipe output",
        "run date",
        "estimated run duration",
        "estimated max memory (GB)",
        "average cpu",
    ]
    lines.append(tr(th(txt) for txt in column_titles))

    for recipe_dir in sorted(Path(output_dir).glob('recipe_*')):
        log = recipe_dir / 'run' / 'main_log.txt'
        success = log.read_text().endswith('Run was successful\n')
        if success:
            status = 'success'
        else:
            debug_log = f"{recipe_dir.name}/run/main_log_debug.txt"
            status = "failed (" + link(debug_log, 'debug') + ")"
        name = recipe_dir.name[:-16]
        date = datetime.datetime.strptime(recipe_dir.name[-15:],
                                          "%Y%m%d_%H%M%S")
        resource_usage = get_resource_usage(recipe_dir)

        entry = []
        entry.append(status)
        entry.append(link(recipe_dir.name, name))
        entry.append(str(date))
        entry.extend(resource_usage)

        entry_txt = tr(td(txt) for txt in entry)
        lines.append(entry_txt)

    return lines


def write_index_html(lines, output_dir):
    """Write lines to index.html."""
    header = textwrap.dedent("""
    <!doctype html>
    <html>
      <head>
        <title>ESMValTool recipes</title>
      </head>
      <style>
      #recipes {
        font-family: Arial, Helvetica, sans-serif;
        border-collapse: collapse;
        width: 100%;
      }

      #recipes td, #recipes th {
        border: 1px solid #ddd;
        padding: 8px;
      }

      #recipes tr:nth-child(even){background-color: #f2f2f2;}

      #recipes tr:hover {background-color: #ddd;}

      #recipes th {
        padding-top: 12px;
        padding-bottom: 12px;
        text-align: left;
        background-color: hsl(200, 50%, 50%);
        color: white;
      }
      </style>
      <body>
        <table id="recipes">
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
