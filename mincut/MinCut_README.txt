## 
## README for NetNC Minimum Cut
##

##@authors Jeremy Owen (jeremyandrewowen@gmail.com) and Ian Overton (i.overton@qub.ac.uk)

## Licensed under the GNU General Public License (GPL) version 3
## you should have recieved a copy of the GNU General Public License
## with the NetNC software (which includes iterative minimum cut)
## in the file LICENSE.txt - if not, see <http://www.gnu.org/licenses/>.

## Installation

Dependencies: numpy, networkx 1.8 (and its dependencies) and optionally matplotlib (for -v option)

All of the required packages should be installable using the command 
line tool 'easy_install', for example, typing 'easy_install networkx==1.8' 
should install the correct version of networkx.

If your system does not find packages installed using easy_install then
you can add their location(s) to the environment variable PYTHONPATH. For example:
setenv PYTHONPATH $PYTHONPATH':'/path/to/python_packages/networkx-1.8-py2.7.egg/

(An example path to python packages is: /usr/lib/python2.7/site-packages/)

## Usage

calling 'python itercut.py' provides the following help message:

USAGE: itercut.py -i [graphfile.txt] -o [outputfile.txt] -t [density threshold value] [-v] [-h]
 Mandatory arguments:
 -i: name of input graph file (table of edges, tab-delimited)
 Optional arguments:
 -o: name of output file (default = input name with _OUT appended)
 -t: threshold density value for iterative cut algorithm (default = 0.65)
 -v: output a PNG showing the graph with nodes colored by predicted pathway (requires matplotlib)
 -c: output a file with the edges cut by the algorithm
 -h: print this help message

Format of input graph files
===============================

Please supply input graphs as text files with each row describing an edge as follows:

nodeID1  nodeID2 weight
nodeID1  nodeID3 weight
nodeID1  nodeID4 weight
nodeID2  nodeID5 weight

where the columns are tab-delimited.


Also, please note that itercut.py occasionally has long runtime (depending on the input graph).
