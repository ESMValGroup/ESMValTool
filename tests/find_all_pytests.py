"""Find all pytests in the tests directory."""
import os


def _find_all_pytests():
    """Return all pytest files in the tests/ dir."""
    pytests_valids = []
    for (dirpath, dirnames, filenames) in os.walk("./tests"):
        ptests = [
            filename for filename in filenames
            if filename.startswith("test_")
            and filename.endswith(".py")]
        for pfile in ptests:
            if os.path.isfile(os.path.join(dirpath, pfile)):
                pytests_valids.append(os.path.join(dirpath, pfile))

    print(pytests_valids)
    return pytests_valids


if __name__ == '__main__':
    _find_all_pytests()
