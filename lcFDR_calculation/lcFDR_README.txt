##
## README for NetNC lcFDR calculation scripts
##

## Licensed under the GNU General Public License (GPL) version 3
## you should have recieved a copy of the GNU General Public License
## with the NetNC software (which includes the local FDR calculation scripts)
## in the file LICENSE.txt - if not, see <http://www.gnu.org/licenses/>.


## @author Ian  Overton (i.overton@qub.ac.uk)


Calculation of node-centric pFDR and summation of Local FDR uses the following scripts:

1. NodeFDR_forNetNC.pl - calculates node-centric pFDR values from the output of NetNC (run the script with no arguments to generate usage information)

2. estimateLocalFDR_conservativeSmoothing.pl - analyses the output of NodeFDR_forNetNC.pl to calculate the sum of local FDR (run the script with no arguments to generate usage information)


