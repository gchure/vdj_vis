// Determine which point mutation is clicked
var mut = endog_sel.value;
endog_filter.group = mut;

// Update the sequence view
seq_view.filters[0] = endog_filter;
seq_source.data.view = seq_view;
seq_source.change.emit();

// Determine point mutant displays
var points = [] // Point mutant identities
var colors = [] // Color for each point mutant
var alphas = [] // Alpha for each point mutant

// Iterate through each point mutant in the sequence source
for (var i = 0; i < seq_source.data['mutant'].length ; i++) {
    // Define the endogenous mutant examined
    var endo = seq_source.data['mutant'][i]
    
    // Determine the point mutant at that position.
    var mutation = seq_source.data['point_mutant'][i]

    // If the endogenous mutation is the same as the selection and is not the 
    // same as the reference, push the point mutant to the array 
    if (endo === mut && mutation !== "No Mutation") {
        points.push(seq_source.data['point_mutant'][i]);
        // Determine if the mutant user is hovering over is the same as the 
        // point mutant or if it is not hovered over at all.
        if (hover_mut === 'None' || hover_mut === '' ||
           hover_mut === 'No Mutation' || hover_mut === undefined ||
           hover_mut === mutation) {        
            // Add the display color to the colors array
            colors.push(seq_source.data['display_color'][i]);

            // Add an alpha of 1 to the alphas array
            alphas.push(1);
           }            

        // If *not* hovered over, make it grey and translucent
        else { 
            colors.push('#c2c2c2');
            alphas.push(0.2);
            }
    }
}
                   


// Update the legend.
var  leg_colors = ['slategrey'];
var  leg_mutant = [mut];
var  leg_alphas = [1];
var  leg_xs = [[-10, -9]];
var  leg_ys = [[-10, -9]];

for (var i = 0; i < colors.length; i++) {
    leg_colors.push(colors[i]);
    if (points[i].includes('-') || points[i].includes('WT')) {
        leg_mutant.push(points[i]);
    }
    else {
    leg_mutant.push(points[i].slice(2))
    }
    leg_alphas.push(alphas[i]);
    leg_xs.push([-10, -9])
    leg_ys.push([-10, -9])
}
leg_source.data['xs'] = leg_xs;
leg_source.data['ys'] = leg_ys;
leg_source.data['c'] = leg_colors;
leg_source.data['mutant'] = leg_mutant;
leg_source.data['alpha'] = leg_alphas;
leg_source.change.emit();

// Update the endogenous percentiles
for (var i = 0; i < endog_percs.length; i++ ) {
    var perc = endog_percs[i];
    endog_percs_view[i].filters[0] = endog_filter;
    perc.data.view = endog_percs_view[i];
    perc.change.emit();
}

// Define arrays of the endogenous plots for easy iteration
var endog_data = [loop_endog,
                  dwell_all_endog, 
                  dwell_cut_endog, 
                  dwell_unloop_endog, 
                  post_endog]
var endog_views = [loop_endog_view,
                   dwell_all_endog_view, 
                   dwell_cut_endog_view, 
                   dwell_unloop_endog_view,
                   post_endog_view]

// Update plots with the endogenous data, views, and endogenous filter
for (var i = 0; i < endog_views.length; i++) { 
   endog_views[i].filters[0] = endog_filter;
   endog_data[i].data.view = endog_views[i]
   endog_data[i].change.emit();}

// Define arrays for iteration over all of the point mutants
var point_data = [loop_point,  
                  dwell_unloop_point,
                  dwell_cut_point, 
                  dwell_all_point];
var point_views = [loop_point_view,
                   dwell_unloop_point_view,
                   dwell_cut_point_view,
                   dwell_all_point_view];
var point_filters = [loop_filter,
                     dwell_unloop_filter, 
                     dwell_cut_filter,
                     dwell_all_filter];
                     

// Update the point percentiles
for (var i = 0; i < point_percs.length; i++ ) {
    var indices = []
    var perc = point_percs[i];
    for (var j = 0; j < point_percs[i].data['mutant'].length; j++) {
        var perc_mut = perc.data['mutant'][j]

        if (points.includes(perc_mut) === true) { 
           indices.push(j) ;
           perc.data['color'][j] = colors[points.indexOf(perc_mut)];
           }
    }
    point_percs_filters[i].indices = indices;
    point_percs_view[i].filters[0] = point_percs_filters[i];
    perc.data.view = endog_percs_view[i];
    perc.change.emit();
}



// Iterate through each point mutant data source
for (var i = 0; i < point_data.length; i++) { 
    // Using index filters for points, so instantiate array to track appropriate
    // index
    var indices = [];

    // Iterate through each entry in the point mutant data source
    for (var j = 0; j < point_data[i].data['mutant'].length; j++) {
        // Define the j-th point mutant
        var point_mutant = point_data[i].data['mutant'][j];

        // If this point mutant is in the `points` array defined from sequence,
        // push the index 
        if (points.includes(point_mutant) === true) { 
            indices.push(j) 

            // Update the displayed color to what is in the color array as well 
            // as the alpha 
            point_data[i].data['color'][j] = colors[points.indexOf(point_mutant)];
            point_data[i].data['alpha'][j] = alphas[points.indexOf(point_mutant)];
        }
     }

    // Update the point views
    point_filters[i].indices = indices;
    point_views[i].filters[0] = point_filters[i];
    point_data[i].data.view = point_views[i];
    point_data[i].change.emit();
}

// To avoid overly complex data structure, populate a blank data source with 
// *only* the data to be displayed as lines. 
var multi_data = [dwell_all_point, dwell_cut_point, dwell_unloop_point, post_point];
var multi_views = [dwell_all_blank, dwell_cut_blank, dwell_unloop_blank, post_blank]
for (var i = 0; i < multi_data.length; i++) { 

    // For multiline, need arrays of arrays
    var xs = [];
    var ys = [];

    // Iterate through each pointmutant that makes up the endogenous mutant
    for (var j = 0; j < points.length; j++) {
        // Instantiate empty arrays for the x and y arrays for the particular 
        // point mutant
        var point_xs = [];
        var point_ys = [];

        // Iterate through each index in the corresponding data source and
        // populate index arrays
        for (var k = 0; k < multi_data[i].data['mutant'].length; k++) { 
            var point_mutant = multi_data[i].data['mutant'][k];
                // If the known point mutant is the same as the point mutant at
                // index k, push the indices, alpha, and color
                if (points[j] === point_mutant) { 
                    point_xs.push(multi_data[i].data['x_val'][k]);
                    point_ys.push(multi_data[i].data['y_val'][k]) ;
            }
               
        }
        // Push the point mutant arrays
        xs.push(point_xs);
        ys.push(point_ys);
    }
    // Update the displayed data source with the correct multiline params.
    multi_views[i].data['xs'] = xs;
    multi_views[i].data['ys'] = ys;
    multi_views[i].data['mutant'] = points;
    multi_views[i].data['c'] = colors;
    multi_views[i].data['alpha'] = alphas;
    multi_views[i].change.emit();
}
