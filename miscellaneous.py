import os
import sys
#import warnings
import logging
import pysam

def ensure_path_replacment(filename,allowReplacment):
	#check that file or folder does not exist
	#if does and allowReplacment warn
	#if does not raise exception
	if os.path.exists(filename):
		if allowReplacment:
			startLog(None)
			logging.warning("File "+filename+" will be replaced")
		else:
			raise Exception("File "+filename+" exists")
	else:
		if not os.path.exists(os.path.dirname(filename)):
			os.makedirs(os.path.dirname(filename))
			
def startLog(filename,loglevel=logging.DEBUG): 
	#TODO
	#Add filename check, so if logger was started with different filename
	#additional handler created and warning shown
	#currently additional call is simply ignored
	
	#if filename != None, logging to console and file
	#otherwise only to console
	
	if not len(logging.getLogger().handlers): 	#check logger was not started
		#start logging
		if filename != None:
			fh=logging.FileHandler(filename=filename) #create log to file handler
			fh.setLevel(logging.NOTSET)
			fh.setFormatter(logging.Formatter('%(levelname)s: %(asctime)s %(message)s')) #format of the file messagaes
			logging.getLogger().addHandler(fh) #bind to main logger

		ch=logging.StreamHandler(sys.stdout) #create log to file handler
		ch.setLevel(logging.NOTSET)
		ch.setFormatter(logging.Formatter('%(levelname)s: %(asctime)s %(message)s')) #format of the file messagaes
		logging.getLogger().addHandler(ch) #bind to main logger
		
		logging.getLogger().setLevel(logging.NOTSET)
		logging.debug("Log started.")
	else:
		if filename != None:
			logging.warning("This call to logging is ignored because there are already associated handlers")
	

# Almost copletely based on http://code.activestate.com/recipes/65287/

def startPmDebug():
	def info(type, value, tb):
	   if hasattr(sys, 'ps1') or not sys.stderr.isatty():
		  # we are in interactive mode or we don't have a tty-like
		  # device, so we call the default hook
		  sys.__excepthook__(type, value, tb)
	   else:
		  import traceback, pdb
		  # we are NOT in interactive mode, print the exception...
		  traceback.print_exception(type, value, tb)
		  print
		  # ...then start the debugger in post-mortem mode.
		  pdb.pm()

	sys.excepthook = info
