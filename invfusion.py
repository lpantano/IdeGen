import sys
import getopt
from optparse import OptionParser
from libs.utils import *
from libs.toolsGenome import createGenome
from libs.alignmentGenome import align
import logging

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv=None):
    usagetxt = "usage: %prog  -p param"
    parser = OptionParser(usage=usagetxt, version="%prog 1.0")
    parser.add_option("-p", "--config", dest="cfile",
                  help="config file", metavar="FILE")
    parser.add_option("-g", "--genome", action="store_true",
                   dest="genome", help="create artificial genome ",default=False)
    parser.add_option("-s", "--star", action="store_true",
                   dest="star", help="align with STAR ",default=False)

    if argv is None:
        argv = sys.argv
    try:
        (options, argv) = parser.parse_args()
        # set up logging to file - see previous section for more details
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            datefmt='%m-%d %H:%M',
                            filename='log',
                            filemode='w')
        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        # set a format which is simpler for console use
        formatter = logging.Formatter('%(asctime)s %(name)-12s: %(levelname)-8s %(message)s')
        # tell the handler to use this format
        console.setFormatter(formatter)
        # add the handler to the root logger
        logging.getLogger('').addHandler(console)

        # Now, define a couple of other loggers which might represent areas in your
        # application:
        con = logging.getLogger('console')
        log = logging.getLogger('file')

        #read config file
        d = readconfig(options.cfile)
        print d
        state = 0
        if options.genome:
            state = createGenome(log,d)
        if options.star and state == 0:
            state = align(d,log)
        # more code, unchanged
        return(state)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, "for help use --help"
        return 2

if __name__ == "__main__":
    sys.exit(main())