import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging

logger = logging.getLogger()


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
    cmap = plt.get_cmap('RdYlGn', max_grade + 2)

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


def do_gcos_table(varname, gcos_expert, gcos_reference):
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

    # Create a list. Each item of the list will be itself a list of strings,
    # corresponding either to the headers or to the GCOS reference values
    contents = list()
    with open(gcos_reference, 'rt') as csvfile:
        s = csv.reader(csvfile, delimiter=",", skipinitialspace=True)
        for row in s:
            contents.append(row)

    # Get the horizontal dimension of the GCOS table
    nx = len(contents[0])

    if not isinstance(gcos_expert, dict):
        assert False, "wrong input type in gcos"
    elif not np.all(np.sort(list(gcos_expert.keys())) == np.sort(contents[0])):
        assert False, "wrong names in gcos"
    else:
        data_contents = []
        for key in contents[0]:
            if gcos_expert[key]["unit"] is None:
                # Read in the product table and store the data.
                # Ignore the header!
                this_unit = ""
            elif gcos_expert[key]["unit"] in ["1", "unkown", "no-unit"]:
                this_unit = ""
            else:
                this_unit = " " + gcos_expert[key]["unit"]
            data_contents.append(str(gcos_expert[key]["value"]) + this_unit)
        contents.append(data_contents)

    ny = len(contents)

    # Check if, possibly, one of the rows of the CSV has not the same number
    # of items
    if sum([1.0 * (len(contents[i]) == nx)
            for i in range(len(contents))]) != ny:
        logger.error("do_gcos_table: unfitting number of columns in " +
                     "reference file")
        raise ValueError("Invalid input: unfitting number of columns in " +
                         "gcos reference file: " + gcos_reference)

    # Create the figure
    fig = plt.figure(figsize=(10, 4))
    # X-Y mesh to plot the color array
    # The Y dimension is written from ny to zero as to write items in visually
    # descending order.
    x_grid, y_grid = np.meshgrid(np.arange(nx + 1), np.arange(ny, -1, -1))

    plt.title("GCOS requirements", fontsize=36)
    # Grid lines
    [plt.plot((0, nx), (y, y), color='k') for y in range(ny)]
    [plt.plot((x, x), (0, ny), color='k') for x in range(nx)]

    for y in range(ny):
        if y == ny - 1:  # if header
            fontweight = "bold"
            fontsize = 20
        else:
            fontweight = "normal"
            fontsize = 20
        for x in range(nx):
            # Read in the "go to line" in the csv and convert it to
            # "go to line" instruction
            # When \n stands in a CSV, python reads \\n
            instring = "\n".join(contents[ny - y - 1][x].split("\\n")).strip()
            plt.text(x + 0.5, y + 0.5, instring, ha='center', va='center',
                     fontweight=fontweight, fontsize=fontsize)

    # Add legend on the left
    plt.text(-0.3, ny - 1.5, "GCOS", rotation=90,
             ha="center", va="center", fontsize=15)
    plt.text(-0.3, ny - 2.5, "ECV\n(averages)", rotation=90,
             ha="center", va="center", fontsize=15)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(0, nx)
    plt.ylim(0, ny)
    plt.tight_layout()

    return fig
