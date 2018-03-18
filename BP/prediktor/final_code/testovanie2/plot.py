import plotly.plotly as py
import plotly.graph_objs as go
import plotly

x = ['Random Forest', 'SVM', 'Ensemble']
y=[0.356,-0.029,0.441]

data = [go.Bar(
x=x,
y=y,
marker=dict(
        color=['rgba(204,204,204,1)', 'rgba(204,204,204,1)', 'rgba(222,45,38,0.8)'])
    )]

layout = go.Layout(
    yaxis=dict(
        title='MCC',
        titlefont=dict(
            size=16,
            color='rgb(107, 107, 107)'
        ),
        tickfont=dict(
            size=14,
            color='rgb(107, 107, 107)'
        )
    ),
    legend=dict(
        x=0,
        y=1.0,
        bgcolor='rgba(255, 255, 255, 0)',
        bordercolor='rgba(255, 255, 255, 0)'
    ),
    barmode='relative',
    bargap=0.15,
    bargroupgap=0.1
)
fig = go.Figure(data=data, layout=layout)
plotly.offline.plot(fig, filename='basic-bar.html')
