#! /usr/bin/env python
import commands,sys,os,subprocess,ROOT,numpy
from optparse import OptionParser
import argparse

eos='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'
basedir='/store/cmst3/user/pharris/'
aparser = argparse.ArgumentParser()
aparser.add_argument('-dir'   ,'--dir'      ,nargs='+',type=str,default='crap')
aparser.add_argument('-base'  ,'--base'     ,nargs='+',type=str,default=['mc_v2'])
aparser.add_argument('-mod'  ,'--mod'       ,nargs='+',type=str,default=[1000])

args = aparser.parse_args()
count=0
basecount=0

def clear():
    global count
    global basecount
    os.system('hadd %s_%s.root *.root ' % (args.dir[0],basecount))
    os.system('mv *.root /tmp/pharris/tmp')
    #os.system('%s cp /tmp/pharris/tmp/%s_%s.root eos/cms/%s/%s/consolidate/' % (eos,args.dir[0],basecount,basedir,args.base[0]))
    os.system('cmsStage /tmp/pharris/tmp/%s_%s.root %s/%s/consolidate/' % (args.dir[0],basecount,basedir,args.base[0]))
    count=0
    basecount=basecount+1

def search(dirname):
    global count
    global basecount
    print 'dir',dirname
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        if label.find('log') > 0 or label == 'failed':
            os.system('%s rm -r %s/%s' % (eos,dirname,label))
            continue
        if label.find('.root') > 0:
            count=count+1
            shortdirname=dirname.replace('eos/cms','')
            os.system('cmsStage %s/%s .' % (shortdirname,label))
            if count % args.mod[0] == 0:
                clear()
            continue
        search('%s/%s' % (dirname,label))

def cleardir(dirname):
    dirSearch = '%s ls %s' %(eos,dirname)
    exists = commands.getoutput(dirSearch)
    for label in exists.splitlines():
        if label.find('.root') == 0:
            continue
        os.system('echo %s rm -r %s/%s >> clearAll.sh' % (eos,dirname,label))
        #os.system('echo %s rm -r %s/%s ' % (eos,dirname,label))
            

os.system('rm -rf /tmp/pharris/tmp/')
os.system('mkdir  /tmp/pharris/tmp/')
os.system('%s mkdir %s/%s/consolidate/' % (eos,basedir,args.base[0]))
search('eos/cms%s/%s/%s' % (basedir,args.base[0],args.dir[0]) )
cleardir('eos/cms%s/%s/%s' % (basedir,args.base[0],args.dir[0]) )
clear()
