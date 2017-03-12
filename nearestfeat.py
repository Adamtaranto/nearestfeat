'''
 ***************************************************************************
 *   (c) Andrew Robinson (andrew.robinson@latrobe.edu.au) 2013             *
 *       La Trobe University &                                             *
 *       Life Sciences Computation Centre (LSCC, part of VLSCI)            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Library General Public License (LGPL)   *
 *   as published by the Free Software Foundation; either version 2 of     *
 *   the License, or (at your option) any later version.                   *
 *                                                                         *
 *   ScienceScripts is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU Library General Public License for more details.                  *
 *                                                                         *
 *   You should have received a copy of the GNU Library General Public     *
 *   License along with ScienceScripts; if not, write to the Free Software *
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
 *   USA                                                                   *
 *                                                                         *
 ***************************************************************************
 
Created on 19/09/2013

@author: arobinson
'''

import sys, getopt

from Bio import SeqIO
from BCBio import GFF



def main(argv):
    '''For each subject feature find the nearest neighbour sequence in each direction.
    
    NOTE: it takes into account the strand of the subject so that its 3' and 5'
    neighbours are reported correctly'''
    
    neighbourGffFilename = None
    subjectGffFilename = None
    subjectFeatureNames = 'gene'
    neighbourFeatureNames = 'region'
    verbose = False
    titles = False
    
    delimiter = "\t"
    naValue = 'N/a'
    
    cmdmsg = '''nearestfeat.py [-h] -S <subject.feature.names> -N <neighbour.feature.names> -s <subject.gff> -n <neighbours.gff>'''
    helpmsg = '''For each subject feature find the nearest neighbour sequence in each direction.
    
%s

 -h        Print this help msg
 -v        Print a summary (counts) at end of processing
 -t        Print titles in first row
 -S <str>  Subject feature names, a comma separated list (no spaces), (default: gene)
 -N <str>  Neighbour feature names, a comma separated list (no spaces), (default: region)
 -s <str>  GFF file containing subject information.
 -n <str>  GFF file containing the Neighbour features.
''' % cmdmsg
    
    # parse arguments
    try:
        opts, args = getopt.getopt(argv,"htvS:N:s:n:",[])
    except getopt.GetoptError:
        print cmdmsg
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print helpmsg
            sys.exit()
        elif opt == "-S":
            subjectFeatureNames = arg
        elif opt == "-N":
            neighbourFeatureNames = arg
        elif opt == "-s":
            subjectGffFilename = arg
        elif opt == "-n":
            neighbourGffFilename = arg
        elif opt == "-v":
            verbose = True
        elif opt == "-t":
            titles = True
        
    # compute extra options
    subjectFeatureNamesList = subjectFeatureNames.split(',')
    neighbourFeatureNamesList = neighbourFeatureNames.split(',')
    
    counter = SafeCounter(['+', '-', '.', '?', 1, -1, 'ne'])
    
    ### (1) parse the neighbour gff file ###
    limitInfo = {'gff_type': neighbourFeatureNamesList}
    neighbourHandle = open(neighbourGffFilename, 'rU')
    neighbourTags = []
    neighbourTagLists = {}
    for sequence in GFF.parse(neighbourHandle, limit_info=limitInfo):
        
        # make a list of neighbour locations
        # Note: this performs gene consolidation via a delayed appending loop
        neighbourTags = []
        lastNeighbour = None
        for feature in sequence.features:
            subject = (int(feature.location.start), int(feature.location.end))
            
            # delayed processing to check for overlaps
            if lastNeighbour: # true on 2+ iteration
                if lastNeighbour[1] >= subject[0]: # overlap?
                    lastNeighbour = (lastNeighbour[0], subject[1])
                else:   
                    neighbourTags.append((lastNeighbour[0], True))
                    neighbourTags.append((lastNeighbour[1], False))
                    lastNeighbour = subject
            else:
                lastNeighbour = subject
            
            counter.inc('ne')
        # ^ next feature ^
        
        # finish off the last neighbour (if exists)
        if lastNeighbour:
            neighbourTags.append((lastNeighbour[0], True))
            neighbourTags.append((lastNeighbour[1], False))
        
        neighbourTags.sort()
        neighbourTagLists[sequence.id] = neighbourTags
        
    # ^ next sequence ^
    neighbourHandle.close()
    
    
    # print column titles if requested
    if titles:
        print str(delimiter).join(['Feature id',
                                   'Label',
                                   "5' Distance",
                                   "3' Distance",
                                   'Strand',
                                  ])
    
    
    ### (2) process the subjects (1-by-1) ###
    limitInfo = {'gff_type': subjectFeatureNamesList}
    subjectHandle = open(subjectGffFilename, 'rU')
    for sequence in GFF.parse(subjectHandle, limit_info=limitInfo):
        
        try:
            neighbourTags = neighbourTagLists[sequence.id]
        except:
            neighbourTags = []
        
        for subject in sequence.features:
            (int(subject.location.start), int(subject.location.end))
            
            i = -1
            
            ## find first tag after subject end ##
            endDist = naValue
            for i in xrange(len(neighbourTags)):
                neighbourTag = neighbourTags[i]
                if neighbourTag[0] >= subject.location.end:
                    if neighbourTag[1]: # is start neighbour?
                        endDist = neighbourTag[0] - subject.location.end
                    else:
                        endDist = 0
                    break
                
            # find the first tag before subject start
            startDist = naValue
            for i in xrange(i,-1,-1): # backwards from where first loop ended
                neighbourTag = neighbourTags[i]
                if neighbourTag[0] <= subject.location.start:
                    if neighbourTag[1]: # is start neighbour?
                        startDist = 0
                    else:
                        startDist = subject.location.start - neighbourTag[0]
                    break
            
            counter.inc(subject.strand)
            label = "%s-%s-%s" % (sequence.id, subject.location.start, subject.location.end)
            fid = label
            if 'ID' in subject.qualifiers:
                fid = " & ".join(subject.qualifiers['ID'])
            if subject.strand in ('+', '.', '?', 1):
                print str(delimiter).join([
                                           fid,
                                           label,
                                           str(startDist),
                                           str(endDist),
                                           "+",
                                           ])
            elif subject.strand in ('-', -1):
                print str(delimiter).join([
                                           fid,
                                           label,
                                           str(endDist),
                                           str(startDist),
                                           "-",
                                           ])
            elif verbose:
                print "Unknown strand: %s" % subject.strand
        
        # ^ next subject ^
    # ^ next sequence ^
    subjectHandle.close()
    
    if verbose:
        print ''
        print '[Feature counts]'
        print 'Total:        %s' % (counter.total - counter.counters['ne'])
        print '+ve:          %s' % (counter.counters['+'] + counter.counters[1])
        print '-ve:          %s' % (counter.counters['-'] + counter.counters[-1])
        print 'non-stranded: %s' % counter.counters['.']
        print 'unknown:      %s' % counter.counters['?']
        print 'others:       %s' % counter.others
        print ''
        print '[Neighbour counts]'
        print 'Total:        %s' % counter.counters['ne']
    
    
class SafeCounter(object):
    ''''''
    
    def __init__(self, lists):
        self.counters = {}
        for val in lists:
            self.counters[val] = 0
        self.total = 0
        self.others = 0
            
    def inc(self, key):
        try:
            self.counters[key] += 1
        except KeyError:
            self.others += 1;
        self.total += 1
## end SafeCounter




if __name__ == "__main__":
    main(sys.argv[1:])

## EOF ##
