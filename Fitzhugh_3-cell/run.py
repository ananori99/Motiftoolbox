#!/usr/bin/env python



#import Tkinter
import system as sys
import network as netw
import traces as tra
import info as nf
import torus as tor
import pylab as pl

#root = Tkinter.Tk()
#screen_width = root.winfo_screenwidth()
#screen_height = root.winfo_screenheight() 

pos_info = '+0+600'
pos_tra = '+300+600'
pos_net = '+300+0'
pos_sys = '+0+0'
pos_torus = '+800+0'

i = nf.info(position=pos_info)
n = netw.network(info=i, position=pos_net)
s = sys.system(info=i, position=pos_sys, network=n)
n.system = s
t = tra.traces(s, n, info=i, position=pos_tra)
tor = tor.torus(s, n, t, info=i, position=pos_torus)

if pl.get_backend() == 'TkAgg':
	s.fig.tight_layout()
	#n.fig.tight_layout()
	t.fig.tight_layout()
	tor.fig.tight_layout()

pl.show()
