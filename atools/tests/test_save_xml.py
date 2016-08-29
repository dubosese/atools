from atools.io import save_hoomdxml
from mbuild.examples.alkane.alkane import Alkane

propane = Alkane(n=3)
box = propane.boundingbox
save_hoomdxml(propane,box,'propane.xml',0,forcefield='opls')
