import streamlit as st
import plotly.express as px
import numpy as np
from src import util

st.set_page_config(page_title = "ClimateDecomposer", page_icon="~")

st.title("ClimateDecomposer")
st.markdown('##### :earth_africa: Decompose paleoclimate records using *Singular Spectrum Analysis* (SSA)  :earth_africa:')

about = '''How to use:
1. Choose a paleoclimate record to decompose from the dropdown menu in the sidebar
2. Choose a **window length** to configure the SSA decomposition
3. Select and view the **Reconstructed Components**: Low RCs (1 to ~10) will usually show the more important components, while high RCs are often noise
4. View the periodogram from the selected RC to check if it aligns with known orbital frequencies.

Created and maintained by [Aidan Starr](https://aidanstarr.github.io/). SSA code adapted from [Braun et al., 2023](https://doi.org/10.1038/s43247-023-00717-5) and [this blog post by Jordon D\'Arcy](https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition/notebook). The app itself is inspired by [**Renewcast** by Giannis Tolios](http://renewcast.giannis.io/) - check it out.
'''

st.markdown(about)

expander = st.expander(":chart_with_upwards_trend: SSA Primer :chart_with_upwards_trend:")
expander.write('''coming soon....
A more comprehensive primer on SSA can be found [in this paper by Vautard and Ghil](https://doi.org/10.1016/0167-2789(89)90077-8).
''')

### sidebar
series = st.sidebar.selectbox(label = "Select a paleo-climate time series",
                               options = ['Benthic Stack','Ice Core CO2'])

container_cat = st.sidebar.container()

with container_cat:
    L = st.sidebar.slider(label = 'Window Size (L)',
                                         min_value = 2, max_value = 600, value = 200)


#### plot the data
st.markdown('---- \n \n \n')
if series == 'Benthic Stack':
    ylabel = 'δ18O (‰)'
    destex = '#### Data: The Global Benthic δ18O stack of [Ahn et al., 2017](https://doi.org/10.1093/climsys/dzx002)'

elif series == 'Ice Core CO2':
    ylabel = 'CO2 (ppm)'
    destex = '#### Data: Composite record of atmospheric CO2 trapped in Antarctic ice cores from [Bereiter et al., 2016](https://doi.org/10.1002/2014GL061957) and references therein'

st.markdown(destex)

df = util.get_data(series)
fig = px.line(x=df.index, y=df.data)
fig.update_traces(line_color='#ef476f')
fig.update_layout(xaxis_title='Time Before Present (kyr))', yaxis_title=ylabel)
st.plotly_chart(fig, use_container_width = True)


#### SSA
ssa = util.SSA(df['data'], L=L)

st.markdown('#### Reconstructed Components from SSA')


with container_cat:
        RC_num0,RC_num1 = st.select_slider(label = 'Reconstructed Component', options = range(1,40),value=[1,1])
        if RC_num0==RC_num1:
            RC_num = RC_num0
        else:
            RC_num = list([RC_num0,RC_num1])

if isinstance(RC_num,int):
    y = ssa.reconstruct(RC_num-1)
    titl = str(RC_num)
else:
    y = ssa.reconstruct(np.arange(RC_num[0]-1,RC_num[1]-1,1))
    titl = '['+str(RC_num0)+':'+str(RC_num1)+']'
fig2 = px.line(x=ssa.reconstruct([0]).index, y=y,title='Reconstructed Component '+titl)
fig2.update_traces(line_color='#ef476f')
fig2.update_layout(xaxis_title='Time Before Present (kyr))', yaxis_title=ylabel)
fig2.update_layout(title=dict(font=dict(size=20), y=1,x=0.6, yref='paper'))
st.plotly_chart(fig2, use_container_width = True)

#### periodogram
st.markdown('#### Periodogram of Reconstructed Components')

detrend = st.checkbox('Detrend?')
if detrend:
    y = y - ssa.reconstruct([0])

f,pxx = util.periodogram(ssa.reconstruct([0]).index,y)


fig3 = px.line(x=f, y=pxx, title='Reconstructed Component '+titl)
fig3.update_traces(line_color='#ef476f')
fig3.update_layout(xaxis_title='Frequency (1/kyr)', yaxis_title='Spectral Power')
fig3.update_xaxes(range=[0, 0.12], row=1, col=1)


fig3.add_vrect(x0=1/120, x1=1/90,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Eccentricity", annotation_position="outside top")
fig3.add_vrect(x0=1/46, x1=1/35,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Obliquity", annotation_position="outside top")
fig3.add_vrect(x0=1/25, x1=1/19,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Precession", annotation_position="outside top")

fig3.update_layout(title=dict(font=dict(size=20), y=0.9,x=0.55, yref='paper'))
st.plotly_chart(fig3, width = 50)


#####
# with st.expander('More plots, please!'):
