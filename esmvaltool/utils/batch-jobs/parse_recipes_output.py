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
    for fil in files:
        with open(fil, "r") as outfile:
            lines = outfile.readlines()
            for line in lines:
                if "Run was successful\n" in line:
                    success_rec.append(fil)

    # typical list elem
    # /home/b/b382109/output_v270/recipe_schlund20jgr_gpp_change_1pct.2378956.out
    recipe_outs = [os.path.basename(ofile) for ofile in success_rec]
    recipe_outs = [f.split(".")[0] + ".yml" for f in recipe_outs]
    
    # return a list of success story recipes
    return recipe_outs


def display_in_md():
    """Print out recipes in Markdown list."""
    todaynow = datetime.datetime.now()
    print(f"## Recipe running session {todaynow}\n")
    with open("all_recipes.txt", "r") as allrecs:
        all_recs = [l.strip() for l in allrecs.readlines()]
    recipe_list = sorted(parse_output_file())
    print("### Succesfully run recipes\n\n")
    print(f"{len(recipe_list)} out of {len(all_recs)} so far\n")
    for rec in recipe_list:
        print("- " + rec)

    # look at fails or still running
    bad_recs = [l for l in all_recs if l not in recipe_list]
    bad_recs = sorted(bad_recs)
    print("\n### Recipes that failed (or are still running)\n")
    print(f"{len(bad_recs)} out of {len(all_recs)} so far\n")
    for rec in bad_recs:
        print("- " + rec)


display_in_md()

