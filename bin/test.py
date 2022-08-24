#External functions
#import modules
import sys, os
import APC as APyC
# import faulthandler
# faulthandler.enable(file=sys.stderr, all_threads=True)
os.environ['KMP_DUPLICATE_LIB_OK']='True'

APyC.Benchmark1000(8)
# faulthandler.dump_traceback(file=sys.stderr, all_threads=True)
print('finished')
