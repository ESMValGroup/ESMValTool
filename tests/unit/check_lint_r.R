library(lintr)
lint("/home/Earth/jvegas/PycharmProjects/ESMValTool/tests/unit/check_lint_r.R")
# check_paths = [
#         'esmvaltool',
#         'tests',
#     ]
#     root_folder = os.path.abspath(os.path.join(__file__, '..', '..', '..'))

#     has_errors = False

#     linters = lintr.with_defaults()

#     for path in check_paths:
#         for dirpath, _, filenames in os.walk(os.path.join(root_folder, path)):
#             for filename in filenames:
#                 if os.path.splitext(filename)[1].lower() == '.r':
#                     errors = lintr.lint(os.path.join(dirpath, filename),
#                                         linters)
#                     for error in errors:
#                         print(str(error)[0:-1])
#                         has_errors = True
    