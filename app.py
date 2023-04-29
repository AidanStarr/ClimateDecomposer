import streamlit as st
import plotly.express as px
from src import util

st.set_page_config(page_title = "ClimateDecomposer", page_icon="~")

st.title("~ClimateDecomposer")
st.markdown('### Decomposing the paleoclimate record into orbital components using Singular Spectrum Analysis')

about = '''ClimateDecomposer allows you to decompose famous paleoclimate records into their constituent components and noise. In other words, it lets you visualise the different components of past climates using *Singular Spectrum Analysis* (SSA), in particular the orbital cycles underlying long-term climate records. A periodogram of each Reconstructed Components can then be viewed to check if it aligns with a known orbital frequency band.  

Created and maintained by [Aidan Starr](https://aidanstarr.github.io/). SSA code adapted from [Braun et al., 2023](https://doi.org/10.1038/s43247-023-00717-5) and [this blog post by Jordon D\'Arcy](https://www.kaggle.com/code/jdarcy/introducing-ssa-for-time-series-decomposition/notebook). The app itself is inspired by [**Renewcast** by Giannis Tolios](http://renewcast.giannis.io/) - check it out.  

A primer on SSA for paleoclimate time series can be found [in this paper by Vautard and Ghil](https://doi.org/10.1016/0167-2789(89)90077-8). For now, try different window sizes (L) using the slider on the left, and see what the different Reconstructed Components are like! 
'''

st.markdown(about)


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


with container_cat:
        RC_num = st.selectbox(label = 'Reconstructed Component', options = range(1,L))   

st.markdown('#### Reconstructed Components from SSA')

fig2 = px.line(x=ssa.reconstruct([0]).index, y=ssa.reconstruct([RC_num-1]),title='Reconstructed Component '+str(RC_num))
fig2.update_traces(line_color='#ef476f')
fig2.update_layout(xaxis_title='Time Before Present (kyr))', yaxis_title=ylabel)
fig2.update_layout(title=dict(font=dict(size=20), y=1,x=0.6, yref='paper'))
st.plotly_chart(fig2, use_container_width = True)

#### periodogram
f,pxx = util.periodogram(ssa.reconstruct([0]).index,ssa.reconstruct([RC_num-1]))

st.markdown('#### Periodogram of Reconstructed Components')

fig3 = px.line(x=f, y=pxx, title='Reconstructed Component '+str(RC_num))
fig3.update_traces(line_color='#ef476f')
fig3.update_layout(xaxis_title='Frequency (1/kyr)', yaxis_title='Spectral Power')
fig3.update_xaxes(range=[0, 0.12], row=1, col=1)


fig3.add_vrect(x0=1/120, x1=1/90,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Eccentricity", annotation_position="outside top")
fig3.add_vrect(x0=1/46, x1=1/35,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Obliquity", annotation_position="outside top")
fig3.add_vrect(x0=1/25, x1=1/19,line_width=0, fillcolor="#dde5b6", opacity=0.1,annotation_text="Precession", annotation_position="outside top")

fig3.update_layout(title=dict(font=dict(size=20), y=1,x=0.6, yref='paper'))
st.plotly_chart(fig3, width = 50)


##### 
# with st.expander('More plots, please!'):
    
