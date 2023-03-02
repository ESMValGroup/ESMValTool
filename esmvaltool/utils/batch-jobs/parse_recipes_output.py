"""Parse recipes output after a mega run."""
import datetime
import os

import glob


def parse_slurm_output(dirname, pattern):
    """Parse the out dir from slurm."""
    pat = os.path.join(dirname, pattern)
    files = glob.glob(pat)

    return files


def parse_output_file():
    """Parse .out files in a given dir."""
    dirname = "/home/b/b382109/output_v270"
    pattern = "*.out*"
    files = parse_slurm_output(dirname, pattern)
    success_rec = []
    diag_fail_rec = []
    missing_data = []
    for fil in files:
        with open(fil, "r", encoding='utf-8') as outfile:
            lines = outfile.readlines()
            for line in lines:
                if "Run was successful\n" in line:
                    success_rec.append(fil)
                elif "esmvalcore._task.DiagnosticError" in line:
                    diag_fail_rec.append(fil)
                elif "ERROR   Missing data for preprocessor" in line:
                    missing_data.append(fil)

    # typical list elem
    # /home/b/b382109/output_v270/recipe_zmnam.2378956.out
    ok_recipe_outs = [os.path.basename(ofile) for ofile in success_rec]
    ok_recipe_outs = [f.split(".")[0] + ".yml" for f in ok_recipe_outs]
    df_recipe_outs = [os.path.basename(ofile) for ofile in diag_fail_rec]
    df_recipe_outs = [f.split(".")[0] + ".yml" for f in df_recipe_outs]
    md_recipe_outs = [os.path.basename(ofile) for ofile in missing_data]
    md_recipe_outs = [f.split(".")[0] + ".yml" for f in md_recipe_outs]

    return (sorted(set(ok_recipe_outs)),
            sorted(set(df_recipe_outs)),
            sorted(set(md_recipe_outs)))


def display_in_md():
    """Print out recipes in Markdown list."""
    todaynow = datetime.datetime.now()
    print(f"## Recipe running session {todaynow}\n")
    with open("all_recipes.txt", "r") as allrecs:
        all_recs = [rec.strip() for rec in allrecs.readlines()]

    # parse different types of recipe outcomes
    recipe_list, failed, missing_dat = parse_output_file()
    print("### Succesfully run recipes\n\n")
    print(f"{len(recipe_list)} out of {len(all_recs)}\n")
    for rec in recipe_list:
        print("- " + rec)

    # surely failed with diagnostic error
    print("\n### Recipes that failed with DiagnosticError\n")
    print(f"{len(failed)} out of {len(all_recs)}\n")
    for rec in failed:
        print("- " + rec)

    # missing data
    print("\n### Recipes that failed of Missing Data\n")
    print(f"{len(missing_dat)} out of {len(all_recs)}\n")
    for rec in missing_dat:
        print("- " + rec)

    # look at other fails or still running
    bad_recs = [
        rec for rec in all_recs
        if rec not in recipe_list and rec not in failed
        and rec not in missing_dat
    ]
    bad_recs = sorted(bad_recs)
    print(
        "\n### Recipes that failed of other reasons  or are still running\n"
    )
    print(f"{len(bad_recs)} out of {len(all_recs)} so far\n")
    for rec in bad_recs:
        print("- " + rec)


if __name__ == '__main__':
    display_in_md()
