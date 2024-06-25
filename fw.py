import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


with st.expander("Viscosities and endpoints", expanded=True):
    col_fluid = st.columns(6)
    with col_fluid[0]:
        muo = st.number_input('$\mu_{o}$', value=1.0, min_value=0.1, max_value=80.0, step=0.1)
    with col_fluid[1]:
        muw = st.number_input('$\mu_{w}$', value=1.0, min_value=0.1, max_value=2.0, step=0.1)
    with col_fluid[2]:
        Krw_sor = st.number_input('$K_{rw}$', value=0.40, step=0.05)
    with col_fluid[3]:
        Kro_swi = st.number_input('$K_{ro}$', value=0.90, step=0.05)
    with col_fluid[4]:
        Swi = st.number_input('$S_{wi}$', value=0.10, step=0.05)
    with col_fluid[5]:
        Sor = st.number_input('$S_{or}$', value=0.22, step=0.05)

c1, c2 = st.columns([2/8, 6/8])
with c1: #corey
    expanderCorey = st.expander("Corey:red_circle:", expanded=True)
    with expanderCorey:
        colCorey = st.columns(2)
        with colCorey[0]:
            no = st.number_input('$n_{o}$', value=2.0, min_value=1.0, max_value=6.0, step=0.1)
        with colCorey[1]:
            nw = st.number_input('$n_{w}$', value=2.0, min_value=1.0, max_value=6.0, step=0.1)

with c2: #let
    expanderLet = st.expander("LET:large_blue_circle:", expanded=True)
    with expanderLet:
        colLet = st.columns(6)
        with colLet[0]:
            lo = st.number_input('$L_{o}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)
        with colLet[1]:
            to = st.number_input('$T_{o}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)
        with colLet[2]:
            eo = st.number_input('$E_{o}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)
        with colLet[3]:
            lw = st.number_input('$L_{w}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)
        with colLet[4]:
            tw = st.number_input('$T_{w}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)
        with colLet[5]:
            ew = st.number_input('$E_{w}$', value=2.0, min_value=1.0, max_value=4.0, step=0.1)

with st.container(border=True):
    col_show = st.columns([1,1,1.5,0.9,1.5,1.5,1.5])
    with col_show[0]:
        show_kr = st.checkbox('$k_{r}$', value=True)
    with col_show[1]:
        show_fw = st.checkbox('$f_{w}$', value=True)
    with col_show[2]:
        show_sec = st.checkbox('secante', value=True)
    with col_show[3]:
        show_der = st.checkbox("f'", value=True)
    with col_show[4]:
        show_Sw = st.checkbox(r'$\bar{S}_{w}$, $S_{wf}$', value=True)
    with col_show[5]:
        show_corey = st.checkbox('Corey:red_circle:', value=True)
    with col_show[6]:
        show_let = st.checkbox('LET:large_blue_circle:', value=True)

Swn = np.linspace(0,1,200)
Sw = Swi + Swn*(1-Sor-Swi) # denorm de sw

# corey;
kro_corey = Kro_swi * np.power((1-Sor-Sw)/(1-Swi-Sor), no)
krw_corey = Krw_sor * np.power((Sw-Swi)/(1-Swi-Sor), nw)
fw_corey = 1/(1+muw/krw_corey * kro_corey/muo)
deriv_corey = np.gradient(fw_corey, Sw)

# let
kro_let = Kro_swi * np.power(1-Swn,lo) / (np.power(1-Swn,lo) + eo*(np.power(Swn, to)))
krw_let = Krw_sor * np.power(Swn,lw) / (np.power(Swn,lw) + ew*(np.power(1-Swn, tw)))
fw_let = 1/(1+muw/krw_let * kro_let/muo)
deriv_let = np.gradient(fw_let, Sw)

fig = plt.figure(figsize=(8,6))
if show_kr:
    if show_corey:
        plt.plot(Sw, krw_corey, 'r', label='Kr Corey', lw=.3)
        plt.plot(Sw, kro_corey, 'r', lw=.3)
    if show_let:
        plt.plot(Sw, krw_let, 'b', label='Kr LET', lw=.3)
        plt.plot(Sw, kro_let, 'b', lw=.3)
if show_fw:
    if show_corey:
        plt.plot(Sw, fw_corey, 'r', label='fw Corey')
    if show_let:
        plt.plot(Sw, fw_let, 'b', label='fw LET')

plt.grid()
plt.ylabel('$K_{ro}, K_{rw}, f_{w}$')
plt.xlabel('$S_{w}$')
plt.xticks(np.arange(0,1.1,0.1))
plt.xlim([0,1])
plt.ylim([0,1.01])
# plt.legend()

s0 = Sw[0]
# corey
for i, s in enumerate(Sw):
    der = deriv_corey[i]
    sec = fw_corey[i]/(s-s0)
    if (sec-der) > 0:
        break
a_corey = fw_corey[i]/(s-s0)
b_corey = -s0*fw_corey[i]/(s-s0)
sm_corey = (1-b_corey)/a_corey

# let
for j, sj in enumerate(Sw):
    der = deriv_let[j]
    sec = fw_let[j]/(sj-s0)
    if (sec-der) > 0:
        break
a_let = fw_let[j]/(sj-s0)
b_let = -s0*fw_let[j]/(sj-s0)
sm_let = (1-b_let)/a_let

# box = dict(pad=2.0, fc=fc, edgecolor='white')
box = dict(boxstyle="round", fc="0.9")

if show_sec:
    if show_corey:
        plt.plot(Sw,a_corey*Sw+b_corey, 'r-.', lw=0.5)
    if show_let:
        plt.plot(Sw,a_let*Sw+b_let, 'b-.', lw=0.5)
if show_Sw:
    if show_corey:
        plt.plot([s], [fw_corey[i]], 'ro')
    if show_let:
        plt.plot([sj], [fw_let[j]], 'bo')
    if show_corey:
        plt.plot([sm_corey], [1], 'ro')
    if show_let:
        plt.plot([sm_let], [1], 'bo')
    if show_corey:
        plt.annotate(f'{s:.2f}', [s-0.03,fw_corey[i]-0.05], c='r')
        plt.annotate(f'{sm_corey:.2f}', [sm_corey-0.03,0.95], c='r')
    if show_let:
        plt.annotate(f'{sj:.2f}', [sj-0.03,fw_let[j]-0.05], c='b')
        plt.annotate(f'{sm_let:.2f}', [sm_let-0.03,0.95], c='b')

if show_der:
    plt.twinx()
    plt.ylabel("$f'_{w}$")
    if show_corey:
        plt.plot(Sw, deriv_corey, 'm')
    if show_let:
        plt.plot(Sw, deriv_let, 'c')
plt.annotate(r'$M^{o}=$'+f'{Krw_sor/Kro_swi*muo/muw:.2f}', [0.02,0.2], bbox=box)
st.pyplot(fig.figure, clear_figure=True)


fig2 = plt.figure(figsize=(8,2))

plt.plot([.5,1], [Swi,Swi], 'b')
plt.plot([.5,.5], [Swi,s],'b')
# plt.plot([0,0.5], [1-Sor,s], 'b')

plt.plot([0,1], [s,s], 'r--')
plt.annotate(r'$S_{wf}$='+f'{s:.2f}', [0.7,s-0.1], c='r')
plt.plot([0,1], [sm_corey,sm_corey], 'k--')
plt.annotate(r'$\bar{S}_{w}$='+f'{sm_corey:.2f}', [0.7,sm_corey+0.05], c='k')

plt.plot([0,1], [1-Sor,1-Sor], 'g--')
plt.annotate(r'1-$S_{or}$='+f'{1-Sor:.2f}', [0.7,1-Sor+0.05], c='g')

plt.plot([0.4,0.4], [sm_corey,1-Sor], 'g', lw=3)

x = np.linspace(0,0.5,20)
k=-2*np.log(s/(1-Sor))
f=(1-Sor)*np.exp(-k*x)
plt.plot(x, f, 'b')

plt.xlim([0,1])
plt.ylim([0,1])
plt.ylabel('$S_{w}$')
plt.xlabel('$x$')
st.pyplot(fig2.figure, clear_figure=True)
