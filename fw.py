import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

c1, c2, c3, c4, c5, c6, c7, c8 = st.columns(8)

with c1:
    muo = st.slider('$\mu_{o}$', value=1.0, min_value=0.1, max_value=2.0, step=0.1)
with c2:
    muw = st.slider('$\mu_{w}$', value=1.0, min_value=0.1, max_value=2.0, step=0.1)
with c3:
    no = st.slider('$n_{o}$', value=1.6, min_value=1.0, max_value=5.0, step=0.2)
with c4:
    nw = st.slider('$n_{w}$', value=1.2, min_value=1.0, max_value=5.0, step=0.2)
with c5:
    Swi = st.number_input('$S_{wi}$', value=0.15)
with c6:
    Sor = st.number_input('$S_{or}$', value=0.30)
with c7:
    Krw_sor = st.number_input('$K_{rw}$', value=0.35)
with c8:
    Kro_swi = st.number_input('$K_{ro}$', value=0.90)

# coef_ang = st.slider(r'$\alpha$', value=0.68, min_value=0.1, max_value=1.0, step=0.01)

Sw = np.linspace(Swi,1-Sor,50)
krw = Krw_sor * np.power((Sw-Swi)/(1-Swi-Sor), nw)
kro = Kro_swi * np.power((1-Sor-Sw)/(1-Swi-Sor), no)

fw = 1/(1+muw/krw*kro/muo)
deriv = np.gradient(fw, Sw)


fig = plt.figure()
plt.plot(Sw, krw, label='krw')
plt.plot(Sw, kro, label='kro')
plt.plot(Sw, fw, 'k', label='fw')
plt.grid()
plt.ylabel('$K_{ro}, K_{rw}, f_{w}$')
plt.xlabel('$S_{w}$')
plt.xticks(np.arange(0,1.1,0.1))
plt.xlim([0,1])
plt.ylim([0,1.01])
plt.legend()

s0 = Sw[0]
for i, s in enumerate(Sw):
    der = deriv[i]
    sec = fw[i]/(s-s0)
    if (sec-der) > 0:
        break

a = fw[i]/(s-s0)
b = -s0*fw[i]/(s-s0)
sm = (1-b)/a
plt.plot(Sw,a*Sw+b, 'k--')
plt.plot([s], [fw[i]], 'ro')
plt.annotate(r'$S_{wf}$='+f'{s:.3f}', [s,fw[i]-0.05])
plt.plot([sm], [1], 'ro')
plt.annotate(r'$\bar{S}_{w}$='+f'{sm:.3f}', [sm-0.16,0.95])

st.pyplot(fig.figure, clear_figure=True)
