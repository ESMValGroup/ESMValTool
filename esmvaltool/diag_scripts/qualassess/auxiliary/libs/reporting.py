import os
import shutil
import subprocess

from .MD_old.METAdata import METAdata


def do_report(report_data, report_title, work_dir, logger,
              signature="SIG", ecv="ECV", dataset="DATASET", latex_opts=False):
    """
    - report_data  a dictionary of a list of plot file names  (.png, ...)
                    including *full path*
                    OR a dictionary containing strings
    - report_title is a string used as title for the report

    Updated: February 18th, 2018 (B. Mueller)
    Updated: July 31st, 2018 (B. Mueller)
    """

    # define output and temporary directories for Sphinx (source, build);
    # add process id to temporary directory names to allow for
    # execution of multiple instances in parallel

    pid = str(os.getpgid(0))

    path_out = work_dir + os.sep + "reporting"    # the final pdf directory
    src_dir = path_out + os.sep + "source_" + \
        report_title.split()[0].lower() + "_" + (
            "_".join(signature) if isinstance(
                signature, list) else signature) + "_" + pid
    # Sphinx source code directory
    bld_dir = path_out + os.sep + "build_" + \
        report_title.split()[0].lower() + "_" + (
            "_".join(signature) if isinstance(
                signature, list) else signature) + "_" + pid
    # Sphinx build directory

    # create output and temporary directories
    if not (os.path.exists(path_out)):
        os.makedirs(path_out)
    if not (os.path.exists(src_dir)):
        os.makedirs(src_dir)
    if not (os.path.exists(bld_dir)):
        os.makedirs(bld_dir)

    # define filename of Sphinx source code file (.rst format)
    output_file = src_dir + os.sep + "report.rst"
    outfile = open(output_file, "w")

    # title (headline)
    my_title = report_title
    outfile.write(my_title + "\n")
    outfile.write("=" * len(my_title) + "\n\n")

    # search for text and plots instances in report_data

    if isinstance(report_data, dict):

        if "listtext" in list(report_data.keys()):

            # if the input data are a dictionary, we create a bullet point list
            # from all key value pairs

            this_title = "Short Information"
            outfile.write(this_title + "\n")
            outfile.write("=" * len(this_title) + "\n\n")

            if isinstance(report_data["listtext"], dict):

                # if the input data are a list of filenames (plots),
                # we simply put all figures with their corresponding caption
                # read from the plot meta data into the report

                for key in report_data["listtext"]:
                    if isinstance(report_data["listtext"][key], dict):
                        outfile.write("* " + key + "\n\n")
                        for key2 in report_data["text"][key]:
                            if isinstance(
                                    report_data["listtext"][key][key2], dict):
                                outfile.write("  * " + key2 + "\n\n")
                                for key3 in report_data["listtext"][key][key2]:
                                    outfile.write(
                                        "    * " + key3 + ": " +
                                        report_data["listtext"]
                                        [key][key2][key3] + "\n")
                                outfile.write("\n")
                            else:
                                outfile.write(
                                    "  * " + key2 + ": " +
                                    report_data["listtext"][key][key2] + "\n")
                        outfile.write("\n")
                    else:
                        outfile.write(
                            "* " + key + ": " +
                            report_data["listtext"][key] + "\n")

                outfile.write(".. raw:: latex \n\n")
                outfile.write("   \clearpage \n")

            else:
                logger.error(
                    "Wrong format in text entry, nothing can be written!")

        else:
            logger.info("No writable list found! " +
                        "There was no 'listtext' in the dictionary!")

        if "freetext" in report_data.keys():

            this_title = "Description"

            if isinstance(report_data["freetext"], (str)):

                if os.path.isfile(report_data["freetext"]):
                    with open(report_data["freetext"], "r") as freetext:
                        with open(os.path.dirname(os.path.abspath(__file__)) +
                                  os.sep + "libs/predef/empty.txt", "r") \
                                as empty:
                            text = freetext.read()
                            if len(set(text) - set(empty.read())):
                                outfile.write(this_title + "\n")
                                outfile.write("=" * len(this_title) + "\n\n")
                                outfile.write(text)
                                outfile.write("\n\n")

                                outfile.write(".. raw:: latex \n\n")
                                outfile.write("   \clearpage \n")

                            else:
                                logger.info("There is still the empty " +
                                            "description from empty.txt!")
                else:
                    outfile.write(this_title + "\n")
                    outfile.write("=" * len(this_title) + "\n\n")
                    outfile.write(report_data["freetext"] + "\n\n")

                    outfile.write(".. raw:: latex \n\n")
                    outfile.write("   \clearpage \n")

            else:
                logger.error(
                    "Wrong format in text entry, nothing can be written!")

        else:
            logger.info("No writable text found! " +
                        "There was no 'freetext' in the dictionary!")

        if "plots" in list(report_data.keys()):

            this_title = "Figure(s)"
            outfile.write(this_title + "\n")
            outfile.write("=" * len(this_title) + "\n\n")

            if isinstance(report_data["plots"], list):
                MD = METAdata()

                # add all plots in plot_list to Sphinx source code file;
                # the figure captions are extracted from the exif file headers
                # of the .png files (if present)

                for f in report_data["plots"]:

                    filename = f.rpartition(os.sep)[-1]
#                    filepath = f.rpartition(os.sep)[0]

                    shutil.copy(f, src_dir)

                    try:
                        caption = MD.read(f).get_dict()[
                            'ESMValTool']['caption']
                    except BaseException:
                        caption = ["Error: empty caption in " + filename]

                    outfile.write(".. figure:: " + filename + "\n"
                                  "   :align:   center" + "\n"
                                  "   :width:   95%" + "\n\n"
                                  "   " + caption[0] + "\n"
                                  )
                    outfile.write(".. raw:: latex \n\n")
                    # does not react on this
                    outfile.write("   \FloatBarrier \n")

            else:
                # TODO?: ERROR function
                # os.path.dirname(os.path.abspath(__file__)) + os.sep +
                logger.error(
                    "Wrong format in plots entry, nothing can be written!")

        else:
            logger.info("No plottable links found! " +
                        "There was no 'plots' in the dictionary!")

        if "freetext" not in list(
                report_data.keys()) and "listtext" not in report_data.keys() \
                and "plots" not in list(report_data.keys()):
            logger.error("Nothing to write reports from!! " +
                         "There was neither 'plots' nor 'text' in " +
                         "the dictionary!")
            return

    outfile.close()
    # copy Sphinx configuration and index file to temporary source directory
    with open(os.path.dirname(os.path.abspath(__file__)) + os.sep +
              "reporting/source/conf.py", 'r') as conf_py:
        content = conf_py.read()
        content = content.replace("++SUBREPORTTITLE++", report_title.title())
        content = content.replace("++VARIABLE++", ecv.title())
        content = content.replace("++DATASET++", dataset)
        with open(src_dir + os.sep + "conf.py", 'w') as target_conf_py:
            target_conf_py.write(content)

    shutil.copy(os.path.dirname(os.path.abspath(__file__)) +
                os.sep + "reporting/source/index.rst", src_dir)

    # set environment variables for Sphinx (source directory and build
    # directory)
    os.environ['SOURCEDIR'] = src_dir
    os.environ['BUILDDIR'] = bld_dir

    if latex_opts is not None:
        # run Sphinx to create a pdf
        oldpath = os.getcwd()
        os.chdir(
            os.path.dirname(
                os.path.abspath(__file__)) +
            os.sep +
            "reporting")
        if not latex_opts:
            with open(os.devnull, 'wb') as devnull:
                subprocess.call("make latexpdf", shell=True,
                                stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.call("make latexpdf", shell=True)

        os.chdir(oldpath)

        # move pdf to the output directory and rename to report_xxx.pdf
        pdfname = path_out + os.sep + "report_" + \
            report_title.split()[0].lower() + "_" + ("_".join(signature) if
                                                     isinstance(signature,
                                                                list)
                                                     else signature) + \
            ".pdf"
        os.rename(
            bld_dir +
            os.sep +
            "latex" +
            os.sep +
            "ESMValToolC3S_511Report.pdf",
            pdfname)

        # clean up temporary directories
        if os.path.exists(src_dir):
            # remove if exists
            shutil.rmtree(src_dir)
            pass
        if os.path.exists(bld_dir):
            # remove if exists
            shutil.rmtree(bld_dir)
            pass

        logger.info("Successfully created " + pdfname + "!")

    else:
        shutil.copy(
            os.path.dirname(
                os.path.abspath(__file__)) +
            os.sep +
            "reporting/PDF_producer.py",
            path_out)
        shutil.copy(
            os.path.dirname(
                os.path.abspath(__file__)) +
            os.sep +
            "reporting/make.bat",
            path_out)
        shutil.copy(
            os.path.dirname(
                os.path.abspath(__file__)) +
            os.sep +
            "reporting/Makefile",
            path_out)
        logger.info("Successfully created PDF_producer.py in " +
                    path_out + "! " +
                    "For report production, please port this directory " +
                    "including all files to a suitable machine.")
