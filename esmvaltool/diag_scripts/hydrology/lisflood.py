"""LISFLOOD diagnostic."""

from esmvaltool.diag_scripts.shared import run_diagnostic

def main(cfg):
    """Process data for use as input to the LISFLOOD hydrological model
    """
    print(cfg)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)