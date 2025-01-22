from pathlib import Path


def file_exists(file_path):
    run_directory = Path.cwd()
    site_recipes_file_path = run_directory / file_path
    return site_recipes_file_path.is_file()
