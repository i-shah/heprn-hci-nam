import os,sys
import subprocess
import shlex
import copy as cp
from functools import reduce

PY_MODS = {}

def copy(x):
    return cp.copy(x)
    
def flatten(L):
    if not L: return []
    if type(L) != list or len(L)==1: return L
    def aaa(a,b):
        if type(a)!=list: a=[a]
        if type(b)!=list: b=[b]
        return a+b

    return [i for i in list(set(reduce(aaa,L))) if i]
    
    
def ifthen(cond,if_true,if_false):
    if cond:
        return if_true
    else:
        return if_false

def shell_command_by_line(cmd):    
    p = subprocess.Popen(shlex.split(cmd),stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while(True):
      retcode = p.poll() #returns None while subprocess is running
      line = p.stdout.readline()
      yield line
      #if(retcode is not None):
      #  break

def run_shell_command(cmd):
    """
    for line in run_shell_command('ls -l'):
        print line
    """
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')    
    

import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart

def send_email(txt="text",subj="Job complete",email='shah.imran@epa.gov'):
    s = smtplib.SMTP('localhost')
    msg = MIMEText(txt)
    H  = dict(Subject=subj,From=email,To=email)
    for k,v in H.items(): msg[k]=v
    s.sendmail(msg['From'], [msg['To']], msg.as_string())
    s.quit()

