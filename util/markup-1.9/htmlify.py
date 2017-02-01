import markup
from markup import oneliner as e
from operator import attrgetter
import exceptions
import os, errno
import pdb
import re
project_type_name = re.compile(".*(CMIP5_fx|OBS|CMIP5)_([^_]+)_.*")



class htmlifyError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class targetError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class figure:
    def __init__(self, diag, fig):
        self.regex_pattern = re.compile(diag + "_([^_]+)_([^_]+)_([^_.]+).*")
        self.diag = diag
        self.fig = fig
        self.var = self.get_var(fig)
        self.field = self.get_field(fig)
        self.descr = self.get_descr(fig)

    def get_var(self, fig):
        var = self.regex_pattern.search(fig).group(1)
        return var

    def get_field(self, fig):
        field = self.regex_pattern.search(fig).group(2)
        return field

    def get_descr(self, fig):
        descr = self.regex_pattern.search(fig).group(3)
        return descr


class Panel(figure):
    def __init__(self, diag, fig_overview):
        figure.__init__(self, diag, fig_overview)


class Single(figure):
    def __init__(self, diag, fig_single):
        figure.__init__(self, diag, fig_single)

        self.name = self.get_name(self.fig)

    def get_name(self, fig):
        name = project_type_name.search(fig).group(2)
        return name

class Files:
    def __init__(self, target):
        self.files = []
        self.target = target

    def append(self, file):
        if (re.search("debug", file.descr) is not None and
            re.search("debug", self.target) is not None):
            self.files.append(file)

        elif ((re.search("debug", file.descr) is None) and
              (re.search("debug", self.target) is None)):
            self.files.append(file)

        else:
            pass

    def sort(self, sort_key):
        self.files.sort(key=attrgetter(sort_key))

class Panels(Files):
    def __init__(self, target):
        Files.__init__(self, target)

class Singles(Files):
    def __init__(self, target):
        Files.__init__(self, target)

class Diags:
    def __init__(self, root_path, diag, panels, singles):
        self.diag = diag
        self.panels = panels
        self.singles = singles
        self.root_path = root_path


def link_to_single(page, singles, panel):
    # Produce link to single plots
    page.p("Individual model plots:", style='padding-top:100px;')
    page.ul()
    if panel.descr == 'JJAS-rp-timelat':
        use_these = singles
    else:
        use_these = [ss for ss in singles if ss.descr == re.sub("(?:-diff)?(?:-page)?[0-9]?", "", panel.descr)]
    
    #remove_rp_timelon = [ss for ss in singles if re.sub("-rp-timelon", "", ss.descr) == panel.descr]
    #diff_removed = [ss for ss in remove_rp_timelon if ss.descr == re.sub("-diff", "", panel.descr)]
#    if panel.descr == 'JJAS-rp-timelat':
#        pdb.set_trace()
    for s in sorted(use_these, key=lambda name:name.name):
        single_html = s.fig + panel.fig + '.html'
        page.li(e.a(s.name, href=single_html))
    page.ul.close()


def write_single(target, non_target, diag, singles, panel, curr_page, run_path, fig_path, id):
    for s in [single for single in singles]:
        single_html = os.path.join(target, s.fig + id + '.html')

        # Create the actual single plot
        curr_page = markup.page(mode='strict_html')
        curr_page.init(css='htmlify.css')
        curr_page.img(src=os.path.join("..", fig_path, diag, s.fig), style='float: left;')

        # Produce links to other single plots
        link_to_single(curr_page, singles, panel)
        curr_page.br()

        for link in non_target:
            rewrite_linked_filename = globals()[target + "_" + link]()
            link_single = rewrite_linked_filename.switch(target, link, s.fig)
            link_panel = rewrite_linked_filename.switch(target, link, panel.fig)
            curr_page.a("Link to " + link + " version",
                        href=os.path.join("..", link, link_single + link_panel + ".html"),
                        class_='link')
            curr_page.br()
        curr_page.a("Back to overview",
                    href='index.html#' + id, class_='link')
        f = open(single_html, "w")
        for line in curr_page.content:
            f.write("%s\n" % line)
        f.close()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

class Filename_switcher(object):
    """ Class to update filenames/links when linking between different folders
    """
    def __init__(self, **kwargs):
        pass

    def switch(self, target, link, fig):
        return fig

class debug_rampup(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = fig
        ret_fig = re.sub("-debug_CMIP5", "_CMIP5", fig)
        ret_fig = re.sub("-debug", "", ret_fig)
        return ret_fig

class debug_normal(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = fig
        ret_fig = re.sub("-debug_CMIP5", "_CMIP5", fig)
        ret_fig = re.sub("-debug", "", ret_fig)
        return ret_fig

class normal_debug(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = re.sub("_CMIP5", "-debug_CMIP5", fig)
        if (re.search("debug", ret_fig) is None):
            ret_fig = re.sub("\.png", "-debug.png", ret_fig)
        return ret_fig

class debug_amip_normal_amip(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = fig
        ret_fig = re.sub("-debug_CMIP5", "_CMIP5", fig)
        ret_fig = re.sub("-debug", "", ret_fig)
        return ret_fig

class normal_amip_debug_amip(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = re.sub("_CMIP5", "-debug_CMIP5", fig)
        if (re.search("debug", ret_fig) is None):
            ret_fig = re.sub("\.png", "-debug.png", ret_fig)
        return ret_fig

class rampup_debug(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = re.sub("_CMIP5", "-debug_CMIP5", fig)
        if (re.search("debug", ret_fig) is None):
            ret_fig = re.sub("\.png", "-debug.png", ret_fig)
        return ret_fig

class rampup_rampdown(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = fig
        ret_fig = re.sub("2010-2090", "2110-2190", fig)
        ret_fig = re.sub("rampup", "rampdown", ret_fig)
        return ret_fig

class rampdown_rampup(Filename_switcher):
    def __init__(self, **kwargs):
        super(Filename_switcher, self).__init__(**kwargs)
        pass

    def switch(self, target, link, fig):
        ret_fig = fig
        ret_fig = re.sub("2110-2190", "2010-2090", fig)
        ret_fig = re.sub("rampdown", "rampup", ret_fig)
        return ret_fig
