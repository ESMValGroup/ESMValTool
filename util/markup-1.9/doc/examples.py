# This code is in the public domain, it comes
# with absolutely no warranty and you can do
# absolutely whatever you want with it.

"""
This module illustrates the use of the markup.py module.

See http://markup.sourceforge.net/
"""

import sys

try:
    import markup
except:
    print __doc__
    sys.exit( 1 )

items = ( "Item one", "Item two", "Item three", "Item four" )
paras = ( "This was a fantastic list.", "And now for something completely different." )
images = ( "thumb1.jpg", "thumb2.jpg", "more.jpg", "more2.jpg" )

page = markup.page( )
page.init( title="My title", 
css=( 'one.css', 'two.css' ), 
header="Something at the top", 
footer="The bitter end." )

page.ul( class_='mylist' )
page.li( items, class_='myitem' )
page.ul.close( )

page.p( paras )
page.img( src=images, width=100, height=80, alt="Thumbnails" )

print page

print '-'*80

title = "Useless Inc."
header = "Some information at the top, perhaps a menu."
footer = "This is the end."
styles = ( 'layout.css', 'alt.css', 'images.css' )

page = markup.page( )
page.init( css=styles, title=title, header=header, footer=footer )
page.br( )

paragraphs = ( "This will be a paragraph.",
"So as this, only slightly longer, but not much.",
"Third absolutely boring paragraph." )

page.p( paragraphs )

page.a( "Click this.", class_='internal', href='index.html' )
page.img( width=60, height=80, alt='Fantastic!', src='fantastic.jpg' )

print page

print '-'*80

images = ( 'egg.jpg', 'spam.jpg', 'eggspam.jpg' )
    
page = markup.page( case='upper' )

page.div( class_='thumbs' )
page.img( width=60, height=80, src=images, class_='thumb' )
page.div.close( )

print page

print '-'*80

titles = ( 'Best featres of M-theory', 'Best bugs in M-theory', 'Branes and brains' )
universities = ( 'Cambridge', 'MIT', 'Amsterdam' )
dates = ( 'January', 'February', 'March' )

myxml = markup.page( mode='xml' )
myxml.init( encoding='ISO-8859-2' )

myxml.cv.open( )
myxml.talk( titles, university=universities, date=dates )
myxml.cv.close( )

print myxml
                                            
print '-'*80

names =     ( 'Alice', 'Bob', 'Eve' )
positions = ( 'encryption', 'encryption', 'eavesdropper' )
locations = ( 'headquarters', 'headquarters', 'unknown' )

myxml = markup.page( mode='xml', onetags=[ 'person', 'location' ], twotags=[ 'company' ] )
myxml.init( )

myxml.company( name='Secret' )
myxml.person( name=names, position=positions, location=locations )
myxml.location( name=( 'headquarters', 'unknown' ), address=( 'here', 'hmmmm' ) )
myxml.company.close( )

print myxml

print '-'*80

page = markup.page( )

page.p( '', id='empty' )
page.form( )
page.input( type='radio', name='sex', value='male', checked=None )
page.input( type='radio', name='sex', value='female' )
page.form.close( )

print page
