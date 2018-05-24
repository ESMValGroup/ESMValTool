'''
Module to create webpage to display extra plots for area assessments

Must have at least one routine with the interface:

publish(area, runs, page)

where

area is an area_assessment.Area object for the given area assessment
runs is a model_run.Runs object for the given area assessment
page is a webpage.Webpage object for the given area assessment
'''

from utility.markup import oneliner
from auto_assess.picture import Picture

IMAGES = [
('South Asian Monsoon Area Metrics' , 'south'),
('South Asian Monsoon Index Metrics' , 'south_indices'),
('South Asian Monsoon Other Metrics' , 'south_other'),
('East Asian Monsoon Area Metrics' , 'east'),
('East Asian Monsoon Index Metrics' , 'east_indices'),
('East Asian Monsoon Other Metrics' , 'east_other')
]


def publish(page):
    '''
    Routine to create webpage to display extra plots for area assessments
    '''

    for (title, fname) in IMAGES:
        page.h2(title)
        pic1 = Picture('metrics_{}.png'.format(fname))
        pic2 = Picture('metrics_alt_{}.png'.format(fname))
        cell1 = oneliner.td(pic1.publish_oneliner(title='', width='700'))
        cell2 = oneliner.td(pic2.publish_oneliner(title='', width='700'))
        row = oneliner.tr(cell1 + cell2)
        page.table(row, width='100%', align='center')
        page.hr()
