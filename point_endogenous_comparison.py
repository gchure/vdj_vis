#%% 
import numpy as np
import pandas as pd
from bokeh.themes import Theme
import bokeh.io
import bokeh.plotting
from bokeh import events
from bokeh.models import (ColumnDataSource, Div, LinearAxis, CustomJS, Text,
                          CDSView, Grid, GroupFilter, Band, Dropdown, HoverTool,
                          IndexFilter, TapTool, ColorBar, Segment, LinearColorMapper,
                          FixedTicker)
from bokeh.layouts import layout, widgetbox
from bokeh.models.widgets import Select
from bokeh.embed import components
import vdj.io
import vdj.stats
bokeh.plotting.output_file('./point_endogenous_comparison.html', mode='inline')

# Load the data sets
loops = pd.read_csv('data/compiled_loop_freq_bs.csv')
dwell_all_data = pd.read_csv('data/compiled_dwell_times.csv')
post_data = pd.read_csv('data/pooled_cutting_probability_posteriors.csv')
pcuts = pd.read_csv('data/pooled_cutting_probability_summary.csv')
loop_data = pd.read_csv('data/compiled_loop_freq_bs.csv')
pcuts = pd.read_csv('data/pooled_cutting_probability.csv')
pcuts.rename(columns={'n_beads': 'n_loops'}, inplace=True)
post_data['color'] = 'slategrey'

# Restrict the posterior distributions.
post_data = post_data[(post_data['hmgb1']==80) & (post_data['salt']=='Mg')]

# Start by trying to figure out the details of picking the point mutants that
# make up the reference
endog_seqs = vdj.io.endogenous_seqs()
reference = endog_seqs['reference']
nt_idx = vdj.io.nucleotide_idx()

# Generate an empty dataframe and a color list from 5' to 3'
mut_df = pd.DataFrame([])
# colors = bokeh.palettes.Category10_10
colors = bokeh.palettes.Category10_10
for key, val in endog_seqs.items():
    color_idx = 0   
    for i, b in enumerate(reference[1]):
        if b != val[1][i]:
            # Find what bases are present in the endogenous that aren't in ref
            base_letter = val[0][i]
            base_idx = nt_idx[val[0][i]] 

            # Define the name. 
            if i < 7:
                loc = 'Hept'
                pos = i + 1
                _pos = pos

            elif (i >= 7)  & (i < 19):
                loc = 'Spac'
                pos = (i + 1) - 7
                _pos = pos + 7

            else:
                loc = 'Non'
                pos = (i + 1) - 19
                _pos = pos + 19
            name = f'12{loc}{reference[0][i]}{pos}{base_letter}'
            display_name = f'Pos.{_pos} {reference[0][i]}\u2192{base_letter}' 
            color = colors[color_idx]
            # Assemble the mutant dataframe
            mut_info = {'position':i, 'base':base_letter, 'base_idx':base_idx,
                        'point_mutant': name, 'mutant':key, 'color':color, 'display_color':color,
                        'display_name': display_name}
            color_idx += 1
        else:
            mut_info = {'position':i, 'base':reference[0][i], 'base_idx':nt_idx[b],
                        'point_mutant': key, 'mutant':key, 
                        'color':'#c2c2c2', 'display_name':f'{key} sequence'}
        if mut_info['point_mutant'] == 'reference':
           mut_info['point_mutant'] = 'WT12rss'
        if mut_info['point_mutant'] == 'V10-95':
            mut_info['point_mutant'] = 'V10-96'
        # Get the loops and cut statistics
        loop_mut = loop_data[(loop_data['mutant']==mut_info['point_mutant']) & (loop_data['salt']=='Mg') & (loop_data['hmgb1']==80)]
        pcut_mut = pcuts[(pcuts['mutant']==mut_info['point_mutant']) & (pcuts['salt']=='Mg') & (pcuts['hmgb1']==80)]
        n_loops = int(loop_mut['n_loops'].unique())
        n_beads = int(loop_mut['n_beads'].unique())
        n_cuts = int(pcut_mut['n_cuts'].unique())
        mut_info['n_loops'] = n_loops 
        mut_info['n_beads'] = n_beads
        mut_info['n_cuts'] = n_cuts
        # if (key != 'WT12rss') & (key != 'reference') & ('cod' not in key.lower()):
        mut_df = mut_df.append(mut_info, ignore_index=True)
seq_source = ColumnDataSource(mut_df)

dfs = []
for g, d in post_data.groupby('mutant'):
    # Determine the class of mutant.
    if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
        parsed_seq = vdj.io.mutation_parser(g) 
        mut_class = 'point'
    else:
        mut_class = 'endogenous'
    d = d.copy()
    d['class'] = mut_class
    d['alpha'] = 1
    dfs.append(d)
post_data = pd.concat(dfs)

# Split into endogenous and point mutations. 
post_endog = post_data[post_data['class']=='endogenous'] 
post_point = post_data[post_data['class']=='point'] 

# Rename the x and y such that I can loop effectively. 
post_point.rename(columns={'probability':'x_val', 'posterior':'y_val'}, 
                    inplace=True)

# Process the looping data into statistics
loops = loops[(loops['hmgb1']==80) & (loops['salt']=='Mg')]
dfs = []
for g, d in loops.groupby('mutant'):
    d = d.copy()
    # Determine the class of mutant. 
    if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
        parsed_seq = vdj.io.mutation_parser(g) 
        mut_class = 'point'
        position = np.where(vdj.io.mutation_parser(g)['seq_idx'] != vdj.io.endogenous_seqs()['reference'][1])[0]

    else:
        mut_class = 'endogenous'
        position = [29]
    if len(position) == 1:
        d['class'] = mut_class
        d['color'] = 'slategrey'
        d['alpha'] = 1
        d['position'] = position[0]
        dfs.append(d)

loops = pd.concat(dfs)

# Process all datasets
# Keep only the HMGB1 = 80 mM and Mg
dwell_all_data = dwell_all_data[(dwell_all_data['hmgb1']==80) & 
                                (dwell_all_data['salt']=='Mg')]

# Separate into cut, unloop
dwell_cut_data = dwell_all_data[dwell_all_data['cut']==1]
dwell_unloop_data = dwell_all_data[dwell_all_data['cut']==0]

# Iterate through the dwell time data  and compute the ECDFS. 
dfs = []
for source in [dwell_all_data, dwell_cut_data, dwell_unloop_data]: 
    bin_dfs = []
    for g, d in source.groupby('mutant'):
        _df = pd.DataFrame()
        # Determine the mutant
        x, y = np.sort(d['dwell_time_min'].values), np.arange(0, len(d), 1) / len(d)
        y[-1] = 1 
        # Stagger the results so I can recreate a step plot
        staircase_y = np.empty(2 * len(d)) 
        staircase_x = np.empty(2 * len(d)) 
        staircase_y[0] = 0
        staircase_y[1::2] = y
        staircase_y[2::2] = y[:-1]
        staircase_x[::2] = x
        staircase_x[1::2] = x

        # Generate another point array so points can be plotted on arrays. 
        point_x = np.zeros(2 * len(d))
        point_y = np.zeros(2 * len(d))
        point_x[1::2] = x
        point_y[1::2] = y
        point_x[point_x == 0] = -1
        point_y[point_y == 0] = -1

        # Assemble the dataframe. 
        _df['x_val'] = staircase_x
        _df['point_x'] = point_x
        _df['y_val'] = staircase_y
        _df['point_y'] = point_y
        _df['mutant'] = g
        _df['color'] = 'slategrey'
        _df['alpha'] = 1
        _df['level'] = 'underlay'
        if ('Spac' in g) | ('Hept' in g) | ('Non' in g):
            _df['class'] = 'point'
        else:
            _df['class'] = 'endogenous'
        bin_dfs.append(_df)     
    dwell_hist = pd.concat(bin_dfs)
    dfs.append(dwell_hist)
dwell_dist, cut_dist, unloop_dist = dfs

# Separate into point  and endogenous dists.  
dwell_all_point = ColumnDataSource(dwell_dist[dwell_dist['class']=='point']) 
dwell_all_endog  = ColumnDataSource(dwell_dist[dwell_dist['class']=='endogenous']) 
dwell_cut_point = ColumnDataSource(cut_dist[cut_dist['class']=='point'])
dwell_cut_endog = ColumnDataSource(cut_dist[cut_dist['class']=='endogenous'])
dwell_unloop_point = ColumnDataSource(unloop_dist[unloop_dist['class']== 'point'])
dwell_unloop_endog = ColumnDataSource(unloop_dist[unloop_dist['class']== 'endogenous'])

# Assemble teh posterior distribution cds
post_endog = ColumnDataSource(post_endog)
post_point = ColumnDataSource(post_point)
post_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
leg_source = ColumnDataSource({'xs': [], 'ys': [], 'c':[], 'mutant':[], 'alpha':[]})

# Make blank dwell time cds for plotting. 
dwell_all_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
dwell_cut_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})
dwell_unloop_blank = ColumnDataSource({'xs':[], 'ys':[], 'c':[], 'mutant':[], 'alpha':[]})

loops['position'] += 1
loop_point = ColumnDataSource(loops[loops['class']=='point'])
loop_endog = ColumnDataSource(loops[loops['class']=='endogenous'])


# Define filters
endog_filter = GroupFilter(column_name='mutant', group='')
target_filter = GroupFilter(column_name='mutant', group='')
dwell_all_filter = IndexFilter(indices=[])
dwell_cut_filter = IndexFilter(indices=[])
dwell_unloop_filter = IndexFilter(indices=[])
loop_filter = IndexFilter(indices=[])
pooled_filter = IndexFilter(indices=[])

# Define the many, many, (many) views
seq_view = CDSView(source=seq_source, filters=[endog_filter])
dwell_all_endog_view = CDSView(source=dwell_all_endog, filters=[endog_filter])
dwell_all_point_view = CDSView(source=dwell_all_point, filters=[dwell_all_filter])

dwell_cut_endog_view = CDSView(source=dwell_cut_endog, filters=[endog_filter])
dwell_cut_point_view = CDSView(source=dwell_cut_point, filters=[dwell_cut_filter])


dwell_unloop_endog_view = CDSView(source=dwell_unloop_endog, filters=[endog_filter])
dwell_unloop_point_view = CDSView(source=dwell_unloop_point, filters=[dwell_unloop_filter])


loop_point_view = CDSView(source=loop_point, filters=[loop_filter])
loop_endog_view = CDSView(source=loop_endog, filters=[endog_filter])

post_endog_view = CDSView(source=post_endog, filters=[endog_filter])
post_point_view = CDSView(source=post_point)

# Define the dropdown for the interactivity
menu_dict = {m:m for m in mut_df['mutant'].unique() if m != 'V10-95'}
menu_dict['DFL161'] = "DFL16.1-5"
menu_dict['DFL1613'] = "DFL16.1-3"
menu = [(v,k) for k, v in menu_dict.items() if k in ['V1-135', 'V8-18', 'V5-43']]
sel = Dropdown(value='', label='Select Endogenous Sequence', menu=menu)

# Define the figure canvas

# Define the figure canvases
seq_ax = bokeh.plotting.figure(width=600, height=50, x_range=[0, 30],
                                y_range=[-0.01, 0.1], tools=[''],
                                toolbar_location=None)
dwell_all_ax = bokeh.plotting.figure(width=275, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], 
                                    x_axis_label='DNA loop lifetime [minutes]',
                                    y_axis_label = 'cumulative distribution', title='all DNA looping events',
                                    y_range=[0, 1],
                                    tools=[''])
dwell_unloop_ax = bokeh.plotting.figure(width=275, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], y_range=[0, 1],
                                    x_axis_label='DNA loop lifetime [minutes]',
                                    y_axis_label = 'cumulative distribution', title='DNA unlooping events',
                                    tools=[''])
                            
dwell_cut_ax = bokeh.plotting.figure(width=275, height=250, x_axis_type='log', 
                                    x_range=[0.5, 80], y_range=[0, 1],
                                    x_axis_label='DNA loop lifetime [minutes]',
                                    y_axis_label = 'cumulative distribution', title='DNA cutting events',
                                    tools=[''])
loop_freq_ax = bokeh.plotting.figure(width=450, height=250,
                                    x_axis_label='initial sequence', 
                                    y_axis_label = 'frequency of DNA loops',
                                    x_range=[0, 32], y_range=[0, 1],
                                    tools=[''], toolbar_location=None)
pcut_ax = bokeh.plotting.figure(width=450, height=250, x_axis_label='probability',
                                y_axis_label='posterior probability',
                                    title = 'probability of DNA cutting',
                                    x_range=[0, 1],
                                tools=[''], toolbar_location=None)

# Define a legend axis and blank it 
leg_ax = bokeh.plotting.figure(width=100, height= 500, tools=[''],
toolbar_location=None, x_range=[0, 1], y_range=[0, 1]) 
bar_ax = bokeh.plotting.figure(width= 350, height=55, tools=[''], 
                            toolbar_location=None)
leg_ax.multi_line('xs', 'ys', color='c', line_width=10, legend='mutant', 
                  source=leg_source, alpha='alpha')
leg_ax.legend.location = 'top_center'
leg_ax.legend.background_fill_color = None
leg_ax.legend.label_text_font_size = '8pt'

# Format the sequence axis to not be colorful.
for ax in [seq_ax, leg_ax, bar_ax]:
    ax.xaxis.visible = False
    ax.yaxis.visible = False
    ax.background_fill_color = None
    ax.outline_line_color = None

# Set the ticker for the x axis o the sequence
ticks = np.arange(1, 30, 1)
ticks[-1] += 2
loop_freq_ax.ray(31, 0, angle=np.pi/2,length=25, line_color='white', line_width=15, alpha=0.75)
loop_freq_ax.xaxis.ticker = ticks
renamed_ticks = {int(t):s for t, s in zip(ticks, list(reference[0]))}
renamed_ticks[31] = 'endo'
loop_freq_ax.xaxis.major_label_overrides = renamed_ticks

# Add A color bar for the confidence intervals
linear_mapper1 = LinearColorMapper(palette=bokeh.palettes.Greys9[1:-2], low=5, high=99)
ticker = FixedTicker(ticks=[10, 25, 50, 75 ,95])
labels = {10:'10%', 25:'25%', 50:'50%', 75:'75%', 95:'95%'}
bar = ColorBar(color_mapper=linear_mapper1, ticker=ticker,
                location=(0, -10), border_line_color=None,
                major_label_overrides=labels, label_standoff=5, width=150, orientation='horizontal',
                height=10, title='confidence interval')
bar_ax.add_layout(bar)

# Add the variant sequences
variant_seq = bokeh.models.glyphs.Text(x='position', y=0, text='base',
                text_color='color', text_font='Courier', text_font_size='26pt')
sequence = seq_ax.add_glyph(seq_source, variant_seq, view=seq_view)

# Assemble the percentiles
endog_percentile_source = []
point_percentile_source = []
point_perc_view = []
endog_perc_view = []
point_perc_filters = []

percs = list(np.sort(loops['percentile'].unique()))
percs.reverse()

loop_freq_ax.triangle(x=31,  y='loops_per_bead', source=loop_endog,
                view=loop_endog_view, fill_color='white', line_color='slategrey', 
                size=6, level='overlay', legend='observed frequency', alpha='alpha')
loop_freq_ax.triangle(x='position', y='loops_per_bead', source=loop_point,
                view=loop_point_view, fill_color='white', line_color='color', 
                size=6, level='overlay', alpha='alpha')
_alphas = [0.2, 0.4, 0.6, 0.8, 0.9, 1]
for i, p in enumerate(percs):
    d_point = loops[(loops['percentile']==p) & (loops['class']=='point')]
    d_endog = loops[(loops['percentile']==p) & (loops['class']=='endogenous')]

    _source_point = ColumnDataSource(d_point)
    _source_endog = ColumnDataSource(d_endog)
    _point_filter = IndexFilter(indices=[]) 
    _view_point = CDSView(source=_source_point, filters=[_point_filter])
    _view_endog = CDSView(source=_source_endog, filters=[endog_filter])
    point_perc_filters.append(_point_filter)
    point_percentile_source.append(_source_point)
    endog_percentile_source.append(_source_endog)
    point_perc_view.append(_view_point)
    endog_perc_view.append(_view_endog)

    point_band = Segment(x0='position', x1='position', y0='low', y1='high', 
                    line_color='color', line_width=10, line_alpha=_alphas[i]) 
    endog_band = Segment(x0=31, x1=31, y0='low', y1='high', 
                    line_color='slategrey', line_width=10, line_alpha=_alphas[i]) 
        
    loop_freq_ax.add_glyph(_source_point, point_band, view=_view_point)
    loop_freq_ax.add_glyph(_source_endog, endog_band, view=_view_endog)

seq_hover_cb = """
var hover_mut = seq_source.data['point_mutant'][cb_data.index['1d'].indices[0]];
"""

sel_cb = """
var hover_mut = 'None';
"""

# Load the callback cod
with open('point_endogenous_comparison.js', 'r') as file:
    unified_code = file.read()
    sel_js = sel_cb + unified_code
    hover_js  = seq_hover_cb + unified_code

js_cbs = []
for cb in [sel_js, hover_js]:
    js_cbs.append(CustomJS(args={'endog_filter':endog_filter, 'leg_source':leg_source,
                    'dwell_all_filter':dwell_all_filter, 'dwell_cut_filter':dwell_cut_filter,
                    'dwell_unloop_filter':dwell_unloop_filter,
                    'loop_filter':loop_filter,
                    'seq_source':seq_source,  'seq_view':seq_view,
                    'dwell_all_endog':dwell_all_endog, 'dwell_all_endog_view': dwell_all_endog_view,
                    'dwell_cut_endog':dwell_cut_endog, 'dwell_cut_endog_view': dwell_cut_endog_view,
                    'dwell_unloop_endog':dwell_unloop_endog, 'dwell_unloop_endog_view': dwell_unloop_endog_view,
                    'dwell_all_point':dwell_all_point, 'dwell_all_point_view': dwell_all_point_view,
                    'dwell_cut_point':dwell_cut_point, 'dwell_cut_point_view': dwell_cut_point_view,
                    'dwell_unloop_point':dwell_unloop_point, 'dwell_unloop_point_view': dwell_unloop_point_view,
                    'dwell_all_blank': dwell_all_blank, 'dwell_cut_blank':dwell_cut_blank,
                    'dwell_unloop_blank':dwell_unloop_blank,
                    'loop_point':loop_point, 'loop_point_view':loop_point_view, 
                    'loop_endog':loop_endog, 'loop_endog_view':loop_endog_view,
                    'post_endog': post_endog, 'post_endog_view':post_endog_view,
                    'post_point':post_point, 'post_view':post_point_view, 
                    'post_endog': post_endog, 'post_endog_view':post_endog_view,
                    'post_blank':post_blank,
                    'endog_sel':sel,
                    'endog_percs': endog_percentile_source,
                    'point_percs': point_percentile_source,
                    'endog_percs_view':endog_perc_view,
                    'point_percs_view': point_perc_view,
                    'point_percs_filters':point_perc_filters,},
        code=cb))

# Plot the ECDFS of the endogenous samples
dwell_all_rend = dwell_all_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_all_endog, view=dwell_all_endog_view,
                   level='overlay')
dwell_cut_rend = dwell_cut_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_cut_endog, view=dwell_cut_endog_view)
dwell_unloop_rend = dwell_unloop_ax.line('x_val', 'y_val', line_width=2, color='slategrey', 
                   source=dwell_unloop_endog, view=dwell_unloop_endog_view)

loop_freq_ax.triangle('position', 'loops_per_bead', source=loop_point, 
                      view=loop_point_view, color='color', alpha='alpha', 
                      size=8)

# Plot the endogenous posterior distributions
pcut_ax.line('probability', 'posterior', source=post_endog, view=post_endog_view,
            line_width=2, color='slategrey', alpha=0.75)

pcut_ax.varea('probability', y1=0, y2='posterior', source=post_endog, view=post_endog_view,
            color='slategrey', alpha=0.25)

# # Plot the dwell distribution for the point mutants. 
dwell_all_ax.multi_line('xs', 'ys', source=dwell_all_blank, color='c',
                        line_width=2, alpha='alpha') 
dwell_cut_ax.multi_line('xs', 'ys', source=dwell_cut_blank, color='c',
                        line_width=2, alpha='alpha') 
dwell_unloop_ax.multi_line('xs', 'ys', source=dwell_unloop_blank, color='c',
                            line_width=2, alpha='alpha')
post_rend = pcut_ax.multi_line('xs', 'ys', source=post_blank, line_width=2, color='c',
                    alpha='alpha')


hover = HoverTool(renderers=[sequence], callback=js_cbs[1],
tooltips=[('mutation', '@display_name'),
          ('number of beads', '@n_beads'),
          ('number of loops', '@n_loops'),
          ('number of cuts', '@n_cuts')])

seq_ax.add_tools(hover)
sel.js_on_change('value', js_cbs[0])

spacer = Div(text='<br/>')

sel_row = bokeh.layouts.row(sel, bar_ax)
sel_col = bokeh.layouts.column(sel_row, seq_ax) 
dwell_row = bokeh.layouts.row(leg_ax, dwell_unloop_ax, dwell_cut_ax, dwell_all_ax)
row1 = bokeh.layouts.row(loop_freq_ax, pcut_ax)
col1 = bokeh.layouts.column(row1, dwell_row)


lay = bokeh.layouts.column(sel_col, col1)

# Set the theme. 
theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#EFEFEF',
                'outline_line_color': '#000000',
            },
            'Axis': {
            'axis_line_color': "black",
            'major_tick_out': 7,
            'major_tick_line_width': 0.75,
            'major_tick_line_color': "black",
            'minor_tick_line_color': "black",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#EFEFEF',
                'border_line_color': '#FFFFFF',
                'border_line_width': 2
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'text_font_style': 'bold',
                'align': 'center',
                'text_font': 'Helvetica',

                'offset': 2,
            }}}

theme = Theme(json=theme_json)
bokeh.io.curdoc().theme = theme

# FOrmat legend details. 
loop_freq_ax.legend.click_policy = 'hide'
bokeh.io.save(lay)



# %%
