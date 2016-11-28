#!/usr/bin/python3
from wand.image import Image
from PyPDF2 import PdfFileReader
from sys import argv

# pdfname= argv[0]
pdfname= 'orbit.pdf'
mypdf= PdfFileReader(pdfname)
pages= mypdf.getNumPages()

