import yaml

from esmvaltool.diag_scripts.shared import run_diagnostic


def main(cfg):
    with open(cfg['setting_name'], 'w') as file:
        yaml.safe_dump(cfg, file)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
