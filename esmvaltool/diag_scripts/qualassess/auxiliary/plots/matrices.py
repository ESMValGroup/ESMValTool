import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging
import re
from collections import OrderedDict

logger = logging.getLogger(os.path.basename(__file__))


def do_smm_table(csv_expert, csv_definitions):
    """
    Authors:  F. Massonnet and B. Hassler
    Creation: February 8th, 2018
    Updated: February 15th, 2018 (B. Mueller)

    Inputs:
        - csv_expert : A CSV file giving the grades assigned to each
                       cell of the System Maturity Matrix (1 to 6)
                       The first line is supposed to be a header with
                       the names of the different categories. This row is not read,
                       but it useful to keep it if the CSV themselves have to be opened
                       by the expert.
                       Fields in the CSV file must be separated by commas.
                       White spaces can be present around the commas.
                       The strings should be quoted (" ... "), especially if there are
                       commas within them.
                       Non relevant information (i.e., where no grade is given) shall be
                       reported by not filling the field.
                       The CSV file should not have an empty line at the end!
                       Example:
                       >> cat input_file.csv
                             "Software Readiness" , "Metadata", "User Documentation", "Uncertainty Characterisation", "Public access feedback and update","Usage"
                              1                   , 4         , 3                   , 2                             , 1                                 , 6
                              3                   , 2         , 1                   , 4                             , 6                                 , 5
                              1                   , 3         , 1                   , 2                             , 4                                 ,
                              6                   ,           , 3                   , 2                             , 6                                 ,

        - csv_definitions : A CSV file (independent of the product to be assessed) giving
                            the standard names for criteria in the System Maturity Matrix.
                            Same conventions as for the csv_expert file apply.
                            The sign "\n" can be used if a linebreak should appear in the
                            name.

                            Example:
                            >> cat definition_file.csv
                               "Software\nReadiness", "Metadata",  "User\nDocumentation", "Uncertainty\nCharacterisation", "Public access,\nfeedback,\nand update", "Usage"
                               "Coding\nStandards"  , "Standards", "Formal description\nof scientific\nmethodology", "Standards", "Public\nAccess/Archive", "Research"
                               "Software\nDocumentation", "Collection\nlevel", "Formal validation\nreport", "Validation", "Version", "Decision\nsupport\nsystem"
                               "Numerical\nReproducibility\nand portability", "File level", "Formal product\nuser guide", "Uncertainty\nquantification", "User\nfeedback",
                               "Security", , "Formal description\nof operations\nconcept", "Automated quality\nmonitoring", "Updates to record" ,


    Output: a .png file including the System Maturity Matrix

    """

    max_grade = 6  # Maximum possible grade.
    # Grades are integers to be given between
    # 1 and max_grade (1, ..., max_grade)

    # Create a list. Each item of the list will be itself a list of strings,
    # corresponding to each word to appear in the System Maturity Matrix
    # (header and words in the matrix itsef)
    definitions = list()
    with open(csv_definitions, 'rt') as csvfile:
        s = csv.reader(csvfile, delimiter=",", skipinitialspace=True)
        for row in s:
            definitions.append(row)

    # Get the dimensions of the System Maturity Matrix, including header.
    # The x-dimension goes horizontally (along columns) while the y-dimension
    # goes vertically (along rows)

    ny = len(definitions)
    nx = len(definitions[0])

    # Check if, possibly, one of the rows of the CSV has not the same
    # number of items
    if sum([1.0 * (len(definitions[i]) == nx)
            for i in range(len(definitions))]) != ny:
        logger.error("do_smm_table: unfitting number of columns in " +
                     "definition file")
        raise ValueError("Invalid input: unfitting number of columns in " +
                         "smm definition file: " + csv_definitions)

    # The grades to be used as color in the System maturity matrix
    grades = np.empty((ny, nx))
    grades[:] = np.nan

    with open(csv_expert, 'rt') as csvfile:
        s = csv.reader(csvfile, delimiter=",", skipinitialspace=True)
        counter_y = 0  # counter to iterate through rows
        for row in s:
            # Check if row length matches definition file
            if counter_y > 0:  # Ignore header
                counter_x = 0
                for item in row:
                    if item.strip() != '':
                        # TO BE IMPROVED WITH NEW BACK END
                        try:
                            grades[counter_y, counter_x] = int(item)
                        except IndexError:
                            logger.error("do_smm_table: " +
                                         "unfitting number of columns in " +
                                         "definition and expert file")
                            raise ValueError("Invalid input: " +
                                             "unfitting number of columns " +
                                             "in smm definition file: " +
                                             csv_definitions +
                                             " and expert file: " + csv_expert)
                    counter_x += 1
            counter_y += 1

    # Create the figure
    fig = plt.figure(figsize=(9, 5))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    # The color map
    cmap = plt.get_cmap('PuBuGn', max_grade + 2)

    # Plot the cells in color
    plt.pcolor(
        x_grid,
        y_grid,
        grades,
        cmap=cmap,
        vmin=-
        0.5,
        vmax=max_grade +
        1.5)
    # Create colorbar at the bottom
    cb = plt.colorbar(
        boundaries=np.arange(
            0.5,
            max_grade + 1),
        ticks=np.arange(
            1,
            nx + 1),
        orientation="horizontal",
        pad=0.05)
    # Remove ticks
    cb.set_ticks([])
    # Put colorbar as the most bottom layer
    cb.ax.zorder = -1

    # Write legend inside colorbar
    [plt.text(0.5 * nx / max_grade + 1.0 * (i - 1) * nx / max_grade,
              -0.65, str(i),
              fontsize=14,
              fontweight="bold",
              ha="center",
              va="center") for i in np.arange(1, max_grade + 1)]

    # Finish polishing the figure
    plt.title("System Maturity Matrix", fontsize=18)
    # Grid lines
    [plt.plot((0, nx), (y, y), color='k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color='k') for x in range(nx)]

    for y in range(ny):
        if y == ny - 1:  # if header
            fontweight = "bold"
            fontsize = 10
        else:
            fontweight = "normal"
            fontsize = 10
        for x in range(nx):
            # Read in the "go to line" in the csv and convert it to
            # "go to line" instruction
            # When \n stands in a CSV, python reads \\n
            instring = "\n".join(
                definitions[ny - y - 1][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha='center', va='center',
                     fontweight=fontweight, fontsize=fontsize)

    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    return fig


def do_gcos_table(varname, gcos_expert, gcos_reference, gcos_expert_file=None):
    """
    Author  :  F. Massonnet
    Creation: February 9th, 2018
    Updated: February 15th, 2018 (B. Mueller)
             March 1st, 2018 (B. Mueller) (no csv for expert)

    Inputs:
        - varname: the variable name to read the right gcos table parts

        - gcos_expert: A CSV file giving the GCOS criteria value for the dataset
                       or product assessed.
                       The first line is supposed to be a header with
                       the names of the different criteria. This row is not read,
                       but it useful to keep it if the CSV themselves have to be opened
                       by the expert.
                       Fields in the CSV file must be separated by commas.
                       White spaces can be present around the commas.
                       The strings should be quoted (" ... "), especially if there are
                       commas within them.
                       Non relevant information (i.e., where no grade is given) shall be
                       reported by not filling the field.
                       The sign "\n" can be used if a linebreak should appear in the
                       name.
                       The CSV file should not have an empty line at the end!

                       Example:
                       >> cat gcos_expert.csv
                            "Accuracy", "Stability", "Frequency",         "Resolution"
                            0.3       , 0.02       , "Hourly\nto weekly", 0.9

        - gcos_reference: A dictionary the reference values (GCOS standards for that product)
                          keys must of course match the expert file!
                          Example:
                          {"Accuracy":{"value":0.3,"unit":"m2"},
                           "Stability":{"value":0.5,"unit":"m2 0.1 y-1"},
                           "Frequency":{"value":28,"unit":"days"},
                           "Resolution":{"value":0.5,"unit":"degrees"}}
    Outputs:
        - A *.png file that includes a two-row table: first row = header,
          second row = numbers from expert with green or red shading depending
          on whether the GCOS standards are met or not

    TODO:
        - is a GCOS standard met if the expert value is less in an absolute sense?
        - how about frequency? reported as "hourly", "decadal", not numbers
    """

    written = False

    try:
        gcos_expert_read = dict()

        with open(gcos_expert_file, 'r') as csv_file:
            reader = csv.reader(csv_file)
            gcos_expert_list = list(reader)

        for it in gcos_expert_list:
            if it[2] in ['']:
                it[2] = None
            else:
                try:
                    it[2] = float(it[2])
                except BaseException:
                    pass
            if it[0] in gcos_expert_read.keys():
                gcos_expert_read[it[0]].update({it[1]: it[2]})
            else:
                gcos_expert_read.update({it[0]: {it[1]: it[2]}})
    except Exception as e:
        logger.error("Calculated GCOS information written at: " +
                     gcos_expert_file +
                     "\nPlease consider corrections, if needed!")
        open(gcos_expert_file, 'a').close()
        with open(gcos_expert_file, 'w') as csv_file:
            writer = csv.writer(csv_file)
            for key, value in gcos_expert.items():
                for key2, value2 in value.items():
                    writer.writerow([key, key2, value2])
        written = True

    if written:
        pass
    elif gcos_expert != gcos_expert_read:
        gcos_expert = gcos_expert_read
        logger.error("Calculated GCOS information overwritten due to " +
                     "differences.\nGCOS information from file: " +
                     gcos_expert_file)
        logger.error("If you think this should not happen, please delete " +
                     "the respective file and run again!")
    else:
        logger.error("No differences found in calculated and read GCOS " +
                     "information from: " + gcos_expert_file)

    contents = list()
    with open(gcos_reference, 'r') as csvfile:
        s = csv.reader(csvfile, delimiter=";", skipinitialspace=True)
        for row in s:
            contents.append(row)

    variable = varname
    #variable = 'Water vapour'

    # convert the list into an array for easier checks for entries
    gcos_ref_array = np.asarray(contents)

    gcos_idx = []
    # look into the third column of the table to find the CMOR ECV name
    for i in range(len(gcos_ref_array[:, 2])):
        if re.search(r'\b' + variable.strip().upper() + r'\b',
                     gcos_ref_array[i, 2].upper()):
            gcos_idx.append(i)

    gcos_contents = OrderedDict()
    gcos_contents["Frequency"] = [""]
    gcos_contents["Resolution"] = [""]
    gcos_contents["Accuracy"] = [""]
    gcos_contents["Stability"] = [""]

    if len(gcos_idx) >= 1:
        gcos_contents["Frequency"] = [gcos_ref_array[gcos_idx[0], 3]]
        gcos_contents["Resolution"] = [gcos_ref_array[gcos_idx[0], 4]]
        gcos_contents["Accuracy"] = [gcos_ref_array[gcos_idx[0], 5]]
        gcos_contents["Stability"] = [gcos_ref_array[gcos_idx[0], 6]]

    # Get the horizontal dimension of the GCOS table
    nx = len(gcos_contents)

    if not isinstance(gcos_expert, dict):
        assert False, "wrong input type in gcos"
    is_equal = np.array_equal(np.sort(list(gcos_expert.keys())),
                              np.sort(list(gcos_contents.keys())))
    if not is_equal:
        assert False, "wrong names in gcos"
    else:
        # data_contents=[]
        for key in gcos_contents:
            if gcos_expert[key]["unit"] is None:
                this_unit = ""
            elif gcos_expert[key]["unit"] in ["1", "unkown", "no-unit"]:
                this_unit = ""
            else:
                this_unit = " " + gcos_expert[key]["unit"]
            gcos_contents[key].append(str(gcos_expert[key]["value"]) +
                                      this_unit)
            ny = len(gcos_contents[key]) + 1

    # Create the figure
    fig = plt.figure(figsize=(12, 7))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    plt.title("GCOS requirements", fontsize=36)
    # Grid lines
    [plt.plot((0, nx), (y, y), color='k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color='k') for x in range(nx)]

    for y in range(ny - 1):
        for (x, key) in enumerate(gcos_contents.keys()):
            instring = "\n".join(gcos_contents[key][y].split("\\n")).strip()
            plt.text(x + 0.5, ny - y - 1.5, instring, ha='center',
                     va='center', fontweight="normal", fontsize=20)
    for (x, key) in enumerate(gcos_contents.keys()):
        instring = "\n".join(key.split("\\n")).strip()
        plt.text(x + 0.5, ny - 0.5, instring, ha='center', va='center',
                 fontweight="bold", fontsize=20)

    # Add legend on the left
    plt.text(-0.3, ny - 1.5, "GCOS", rotation=90,
             ha="center", va="center", fontsize=15)
    plt.text(-0.3, ny - 2.5, "ECV\n(averages)",
             rotation=90, ha="center", va="center", fontsize=15)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    return fig


def do_eval_table(varname, eval_expert, eval_data, ecv_length):
    """
    Author  : Arun Rana
    Creation: July 26th, 2018
    Updated : August 08th, 2018

    Inputs:
        - varname: the variable name to read the right evaluation/validation table parts

        - eval_expert: A CSV file giving the Evaluation/validation criteria value for the dataset
                       or product assessed.
                       The first line is supposed to be a header with
                       the names of the different criteria. This provides
                       the criteria for.
                       Fields in the CSV file must be separated by commas.
                       White spaces can be present around the commas.
                       The strings should be quoted (" ... "), especially if there are
                       commas within them.
                       Non relevant information (i.e., where no validation is given) shall be
                       reported by not filling the field.
                       The sign "\n" can be used if a linebreak should appear in the
                       name.
                       The CSV file should not have an empty line at the end!

                       Example:
                       >> cat eval_expert.csv
                            "ECV", "Mean\n climatology", "Trend", "Variability"
                             tas,	20,	42,	0.2

        - eval_data: A dictionary of values (read from in the product). If the expert input
                          is smaller than the  length of the dataset for ECV, that cell would get the number for "green".
                          If the expert input is bigger than the  length of the dataset for ECV,
                          that cell would get the number for "red".
                          It is a constant value refering to the total length of data that is read in and is pointed in overview report
                          and saved for GCOS requirement as well. (tim_range_spec)
    Outputs:
        - A *.png file that includes a one-row table with green or red or white shading depending
          on whether the evaluation standards are met or not

    TODO:

    """

#    # WILL HAVE TO BE DELETED
#    # -----------------------
#    path_out = work_dir + "/plot_scripts/python/system_maturity_matrix"
#
#    path_report_out = work_dir + "/reports/sphinx/source"

    # ----------------------

    max_grade_eval = 3  # Maximum possible grade. Grades are integers to be given between
    # 1 and max_grade_eval (1, ..., max_grade_eval)

    # Create a list. Items of the list will be itself a list of strings, corresponding to each
    # word to appear in the Evaluation Matrix (header of the matrix itself -
    # parameter for evaluation)
    definitions = list()
    with open(eval_data, 'r') as csvfile:
        s = csv.reader(csvfile, delimiter=",", skipinitialspace=True)
        for row in s:
            definitions.append(row)

    # Get the dimensions of the Validation Matrix, which is just header in this case.
    # The x-dimension goes horizontally (along columns) while the y-dimension goes
    # vertically (along rows)
    ny = len(definitions)
    nx = len(definitions[0])

    # The grades to be used as color in the Validation matrix
    grades = np.empty((ny, nx))
    grades[:] = np.nan

    # perform the comparison with validation data and prepare the grading matrix for plotting
    # here the grades matrix is to be filled with number assigned with colors
    # below it give it right feel

    time_range = list()
    with open(eval_expert, 'r') as csvfile:
        s = csv.reader(csvfile, delimiter=",", skipinitialspace=True)
        for row in s:
            time_range.append(row)

    # excluding the header for comparison with range (was [1][1:])
    time_range = time_range[1][:]

    # making sure that the conversion to integers does not crash the namelist
    try:
        # conversion to integer for comparison to ECV data length
        time_range = [int(i) for i in time_range]
    except BaseException:
        time_range = [np.nan for i in time_range]
    # print time_range - done during coding

    for i in range(len(time_range)):
        # it's ok if the data are just as long as the limit given by the expert
        # user
        if time_range[i] <= ecv_length:
            grades[0][i] = 3
        elif time_range[i] > ecv_length:
            grades[0][i] = 1
        else:
            grades[0][i] = 2
        # print grades.item (i)

    # Create the figure
    fig = plt.figure(figsize=(6, 2))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    # The color map
    cmap = plt.get_cmap('RdYlGn', max_grade_eval)

    # Plot the cells in color
    plt.pcolor(x_grid, y_grid, grades, cmap=cmap, vmin=-0.5,
               vmax=max_grade_eval + 1.5)

    # Finish polishing the figure
    plt.title("ESM Evaluation Table and Grading", fontsize=18)
    # Grid lines
    [plt.plot((0, nx), (y, y), color='k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color='k') for x in range(nx)]

    for y in range(ny):
        if y == ny - 1:  # if header
            fontweight = "bold"
            fontsize = 15
        else:
            fontweight = "normal"
            fontsize = 15
        for x in range(nx):
            # Read in the "go to line" in the csv and convert it to "go to line" instruction
            # When \n stands in a CSV, python reads \\n
            num = ny - y - 1
            instring = "\n".join(definitions[num][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha='center', va='center',
                     fontweight=fontweight, fontsize=fontsize)

    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    return fig
