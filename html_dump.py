__author__ = 'jpresern'

import sys
import os
# sys.path.append("../../../python_projects/markup")

import markup as mup
# from markup import oneliner as one
# from os.path import join
from posixpath import join as ppjoin
from utils import folder_walk
import fnmatch

"""
Looks inside the experimental subfolder for the specified *.svg images
: param path_to_folder : path to folder containing experimental subfolders
: param subfolder : subfolder containing experimental files
: param type : image type we are interested in (FI_curves, VI_curves,...)
: return matches : list of image file matching criteria

"""
def browse_for_images (path_to_folder, subfolder, type):
    matches = []
    (_,_,filenames) = folder_walk(ppjoin(path_to_folder, subfolder))
    for filename in fnmatch.filter(filenames, ".".join([type,'svg'])):
        matches.append(filename)

    return matches

"""
Generates .html page for each experiment, embedding the figures and returning the path to the hmtml
: param path_to_folder : path to folder containing experimental subfolders
: param subfolder : subfolder containing experimental files
: param info : dictionary containing experiment meta from info.dat
: param figs : ordered dictionary containing figure names
: returns html_name : relative link to the generated html page
: returns rec_location : ELL segment in which the recording was made as poached from info.dat
"""

def generate_exp_html (path_to_folder, subfolder, info, figs):

    #   generate page title
    rec_location = info["Cell"]["Location"]
    pagetitle = "".join([rec_location,": ", subfolder])

    #   creates page instance
    pg = mup.page( )

    #   initialize the page, embed css, titles, header and footer
    pg.init( title=pagetitle,
                css=( 'markup.css', 'two.css' ),
                header="".join(['<font size="10" color="red" >',pagetitle,'</font>']),
                footer="The bitter end." )

    #   line break under header
    pg.br( )

    #   loop over the figs ordered dictionary for every type of the figure
    for k, v in figs.items():
        images = []
        pg.p("".join(['<font size="5" color="red" >','<b>',k,'</b>','</font>']))
        for i in range(0, len(v)):
            # images.append (".".join([v[i],'svg']))
            pg.p()
            v[i].split("_FI")[0]
            pg.a(href = "".join([v[i].split("_FI")[0], "_FI_rug.html"]))
            pg.img( src=".".join([v[i],'svg']), alt= "click for rug and return plot", align="middle")
        pg.br()
        pg.hr()

    #   line break before footer
    pg.br()

    #   dump the generated html code into the .html file
    html_name = ".".join(["".join([info["Cell"]["Location"],"_","index"]),"html"])
    t = []
    f = open(ppjoin(path_to_folder, subfolder, html_name), "w")
    t.append(str(pg))
    f.write(t[0])
    f.close()

    return html_name, rec_location

    # print (pg)
"""
Generate .html containing figures of the FI rug & return plots
: param path_to_folder : folder containing the experiments
: param subfolder : subfolder, containing cell measurements
: param info : a dictionary containing meta from experiment info.dat
: param figs : ordered dictionary containing figure names
: return html_name : relative link to the FI rug & return plot subpage

"""
def generate_rug_html (path_to_folder, subfolder, info, fig_rugs):

    #   loop over the figs ordered dictionary for every type of the figure
    for k, v in fig_rugs.items():

        for i in range(0, len(v)):

            #   generate page title
            rec_location = info["Cell"]["Location"]
            pagetitle = "".join([rec_location,": ", subfolder, ", Rug & Return map plot"])

            #   creates page instance
            pg = mup.page( )

            #   initialize the page, embed css, titles, header and footer
            pg.init( title=pagetitle,
                        css=( 'markup.css', 'two.css' ),
                        header="".join(['<font size="10" color="red" >',pagetitle,'</font>']),
                        footer="The bitter end." )

            #   line break under header
            pg.br()

            pg.p("".join(['<font size="5" color="red" >','<b>',k,'</b>','</font>']))

            pg.p()
            pg.img( src=".".join([v[i],'svg']), align="middle")
            pg.br()

            #   dump the generated html code into the .html file
            html_name = ".".join([v[i],"html"])
            t = []
            f = open(ppjoin(path_to_folder, subfolder, html_name), "w")
            t.append(str(pg))
            f.write(t[0])
            f.close()


    # return html_name

"""
Generates index.html from the links in the dictionary
: param exp_pages : dictionary containing links to the pages of a certain experiment

"""

def generate_index (exp_pages, wd, fish_n_chips):

    pagetitle = "eFISH intracellular in vitro"
    pg = mup.page ()
    pg.init ( title=pagetitle,
                css=( 'markup.css', 'two.css' ),
                header="".join(['<font size="10" color="red" >',pagetitle,'</font>']),
                footer="The void." )

    for k, v in exp_pages.items():
        pages = []
        pg.p("".join(['<font size="5" color="red" >','<b>',k,'</b>','</font>']))
        pg.hr()
        pg.br()
        for i in range(0, len(v)):
            if fish_n_chips[v[i].split("/")[-2]] == []:
                # pg.p(fish_n_chips[v[i]])
                continue
            else:
                pg.p(fish_n_chips[v[i].split("/")[-2]][0])
                pg.a("".join([v[i], "     ", str(fish_n_chips[v[i].split("/")[-2]][1:])]), class_='internal', href=v[i] )

        pg.br()
        pg.hr()


    #   dump the generated html code into the .html file
    html_name = "index.html"
    t = []
    f = open(ppjoin(wd, html_name), "w")
    t.append(str(pg))
    # print(t[0])
    f.write(t[0])
    f.close()

    # return html_name, rec_location
