"""Barebone example of executing an interactive diagnostic."""
from esmvaltool.diag_scripts.shared import run_diagnostic_interactive
from IPython import embed


def main(cfg):
    """Function to jump into ipython."""
    embed()


if __name__ == '__main__':

    with run_diagnostic_interactive() as config:
        main(config)
