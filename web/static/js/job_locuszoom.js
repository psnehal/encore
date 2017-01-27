/* eslint-env jquery */
/* eslint no-unused-vars: ["error", { "vars": "local" }] */
/* global LocusZoom, job_id, api_url */

var lzplot;
$(document).ready(function() {

    var data_sources = new LocusZoom.DataSources();
    var EpactsDS = LocusZoom.Data.Source.extend();
    EpactsDS.prototype.getURL = function (state)
    {
        return api_url + "?chrom=" + state.chr + "&start_pos=" + state.start + "&end_pos=" + state.end;
    };
    var EpactsLD = LocusZoom.Data.Source.extend(function(init) {
        this.parseInit(init);
    }, "LDEP", LocusZoom.Data.LDSource);
    EpactsLD.prototype.findMergeFields = function() {
        return {
            id: "epacts:MARKER_ID",
            position: "epacts:BEGIN",
            pvalue: "epacts:PVALUE|neglog10"
        };
    };

    var apiBase = "/api/lz/";
    data_sources.add("epacts", new EpactsDS)
      .add("ld", ["LDEP", apiBase + "ld-"])
      .add("gene", ["GeneLZ", { url: apiBase + "gene", params: {source: 2} }])
      .add("recomb", ["RecombLZ", { url: apiBase + "recomb", params: {source: 15} }])
      .add("constraint", ["GeneConstraintLZ", { url: apiBase + "constraint" }])
      .add("sig", ["StaticJSON", [{ "x": 0, "y": 4.522 }, { "x": 2881033286, "y": 4.522 }] ]);

    LocusZoom.TransformationFunctions.set("scinotation", function(x) {
        var log;
        if (x=="NA") {return "-";}
        x = parseFloat(x);
        if (Math.abs(x) > 1){
            log = Math.ceil(Math.log(x) / Math.LN10);
        } else {
            log = Math.floor(Math.log(x) / Math.LN10);
        }
        if (Math.abs(log) <= 3){
            return x.toFixed(3);
        } else {
            return x.toExponential(2).replace("+", "").replace("e", " × 10^");
        }
    });


    function getlayout(avail_fields) {
        var has_maf = avail_fields.indexOf("MAF") !== -1;
        var has_beta = avail_fields.indexOf("BETA") !== -1;
        var fields = ["epacts:MARKER_ID", "epacts:CHROM", 
            "epacts:END", "epacts:BEGIN", "epacts:PVALUE|neglog10", 
            "epacts:PVALUE|scinotation", "epacts:PVALUE", 
            "epacts:NS", "ld:state", "ld:isrefvar"];
        var tooltip = "<div style='text-align: right'>"
            + "<strong>{{epacts:MARKER_ID}}</strong><br>"
            + "Chrom: <strong>{{epacts:CHROM}}</strong><br/>"
            + "Pos: <strong>{{epacts:BEGIN}}</strong><br/>"
            + "P Value: <strong>{{epacts:PVALUE|scinotation}}</strong><br>"
            + ((has_maf)? "MAF: <strong>{{epacts:MAF}}</strong><br/>" : "")
            + ((has_beta) ? "BETA: <strong>{{epacts:BETA}}</strong><br/>" : "")
            + "N: <strong>{{epacts:NS}}</strong><br/>"
            + "</div>";
        if (has_maf) {
            fields.push("epacts:MAF");
        }
        if (has_beta) {
            fields.push("epacts:BETA");
        }
        return {
            "state": {},
            "width": 1000,
            "height": 500,
            "resizable": "responsive",
            "panel_boundaries": false,
            "aspect_ratio": 1.7777777777777777,
            "min_region_scale": 20000,
            "max_region_scale": 4000000,
            "panels": [{
                "id": "positions",
                "title": "",
                "width": 800,
                "height": 225,
                "origin": {
                    "x": 0,
                    "y": 0
                },
                "min_width": 400,
                "min_height": 200,
                "proportional_width": 1,
                "proportional_height": 0.5,
                "proportional_origin": {
                    "x": 0,
                    "y": 0
                },
                "margin": {
                    "top": 35,
                    "right": 50,
                    "bottom": 40,
                    "left": 50
                },
                "inner_border": "rgba(210, 210, 210, 0.85)",
                "axes": {
                    "x": {
                        "label_function": "chromosome",
                        "label_offset": 32,
                        "tick_format": "region",
                        "extent": "state"
                    },
                    "y1": {
                        "label": "-log\u2081\u2080 p-value",
                        "label_offset": 28
                    },
                    "y2": {
                        "label": "Recombination Rate (cM/Mb)",
                        "label_offset": 40
                    }
                },
                "interaction": {
                    "drag_background_to_pan": true,
                    "drag_x_ticks_to_scale": true,
                    "drag_y1_ticks_to_scale": true,
                    "drag_y2_ticks_to_scale": false,
                    "scroll_to_zoom": false,
                    "x_linked": true
                },
                "data_layers": [{
                    "id": "recomb",
                    "type": "line",
                    "fields": ["recomb:position", "recomb:recomb_rate"],
                    "z_index": 1,
                    "style": {
                        "stroke": "#0000FF",
                        "stroke-width": "1.5px"
                    },
                    "x_axis": {
                        "field": "recomb:position"
                    },
                    "y_axis": {
                        "axis": 2,
                        "field": "recomb:recomb_rate",
                        "floor": 0,
                        "ceiling": 100
                    },
                    "transition": {
                        "duration": 200
                    }
                }, {
                    "id": "positions",
                    "type": "scatter",
                    "point_shape": "circle",
                    "point_size": {
                        "scale_function": "if",
                        "field": "ld:isrefvar",
                        "parameters": {
                            "field_value": 1,
                            "then": 80,
                            "else": 40
                        }
                    },
                    "color": [{
                        "scale_function": "if",
                        "field": "ld:isrefvar",
                        "parameters": {
                            "field_value": 1,
                            "then": "#9632b8"
                        }
                    }, {
                        "scale_function": "numerical_bin",
                        "field": "ld:state",
                        "parameters": {
                            "breaks": [0, 0.2, 0.4, 0.6, 0.8],
                            "values": ["#357ebd", "#46b8da", "#5cb85c", "#eea236", "#d43f3a"]
                        }
                    }, "#B8B8B8"],
                    "fields": fields,
                    "id_field": "epacts:MARKER_ID",
                    "z_index": 2,
                    "x_axis": {
                        "field": "epacts:BEGIN"
                    },
                    "y_axis": {
                        "axis": 1,
                        "field": "epacts:PVALUE|neglog10",
                        "floor": 0,
                        "upper_buffer": 0.1,
                        "min_extent": [0, 10]
                    },
                    "transition": {
                        "duration": 200
                    },
                    "highlighted": {
                        "onmouseover": "on",
                        "onmouseout": "off"
                    },
                    "selected": {
                        "onclick": "toggle_exclusive",
                        "onshiftclick": "toggle"
                    },
                    "tooltip": {
                        "closable": false,
                        "show": {
                            "or": ["highlighted"]
                        },
                        "hide": {
                            "and": ["unhighlighted"]
                        },
                        "html":  tooltip
                    }
                }]
            }, {
                "id": "genes",
                "width": 800,
                "height": 225,
                "origin": {
                    "x": 0,
                    "y": 225
                },
                "min_width": 400,
                "min_height": 112.5,
                "proportional_width": 1,
                "proportional_height": 0.5,
                "proportional_origin": {
                    "x": 0,
                    "y": 0.5
                },
                "margin": {
                    "top": 20,
                    "right": 50,
                    "bottom": 20,
                    "left": 50
                },
                "axes": {},
                "interaction": {
                    "drag_background_to_pan": true,
                    "scroll_to_zoom": true,
                    "x_linked": true
                },
                "data_layers": [{
                    "id": "genes",
                    "type": "genes",
                    "fields": ["gene:gene", "constraint:constraint"],
                    "id_field": "gene_id",
                    "highlighted": {
                        "onmouseover": "on",
                        "onmouseout": "off"
                    },
                    "selected": {
                        "onclick": "toggle_exclusive",
                        "onshiftclick": "toggle"
                    },
                    "transition": {
                        "duration": 200
                    },
                    "tooltip": {
                        "closable": true,
                        "show": {
                            "or": ["highlighted", "selected"]
                        },
                        "hide": {
                            "and": ["unhighlighted", "unselected"]
                        },
                        "html": "<h4><strong><i>{{gene_name}}</i></strong></h4><div style=\"float: left;\">Gene ID: <strong>{{gene_id}}</strong></div><div style=\"float: right;\">Transcript ID: <strong>{{transcript_id}}</strong></div><div style=\"clear: both;\"></div><table><tr><th>Constraint</th><th>Expected variants</th><th>Observed variants</th><th>Const. Metric</th></tr><tr><td>Synonymous</td><td>{{exp_syn}}</td><td>{{n_syn}}</td><td>z = {{syn_z}}</td></tr><tr><td>Missense</td><td>{{exp_mis}}</td><td>{{n_mis}}</td><td>z = {{mis_z}}</td></tr><tr><td>LoF</td><td>{{exp_lof}}</td><td>{{n_lof}}</td><td>pLI = {{pLI}}</td></tr></table><div style=\"width: 100%; text-align: right;\"><a href=\"http://exac.broadinstitute.org/gene/{{gene_id}}\" target=\"_new\">More data on ExAC</a></div>"
                    }
                }]
            }]
        };
    }

    function move(plot, direction) {
        // 1 means right, -1 means left.
        var start = plot.state.start;
        var end = plot.state.end;
        var shift = Math.floor((end - start) / 2) * direction;
        plot.applyState({
            chr: plot.state.chr,
            start: start + shift,
            end: end + shift
        });
    }

    function zoom(plot, growth_factor){
        // 2 means bigger view, 0.5 means zoom in.
        growth_factor = parseFloat(growth_factor);
        var delta = (plot.state.end - plot.state.start) * (growth_factor - 1) / 2;
        var new_start = Math.max(Math.round(plot.state.start - delta), 1);
        var new_end   = Math.round(plot.state.end + delta);
        if (new_start == new_end){ new_end++; }
        var new_state = {
            start: new_start,
            end: new_end
        };
        if (new_state.end - new_state.start > plot.layout.max_region_scale){
            delta = Math.round(((new_state.end - new_state.start) - plot.layout.max_region_scale) / 2);
            new_state.start += delta;
            new_state.end -= delta;
        }
        if (new_state.end - new_state.start < plot.layout.min_region_scale){
            delta = Math.round((plot.layout.min_region_scale - (new_state.end - new_state.start)) / 2);
            new_state.start -= delta;
            new_state.end += delta;
        }
        plot.applyState(new_state);
    }

    LocusZoom.createCORSPromise("GET",api_url).then(function(x) {
        x = JSON.parse(x);
        var cols = [];
        if (x.header && x.header.variant_columns) {
            cols = x.header.variant_columns;
        }
        var layout = getlayout(cols);
        lzplot = LocusZoom.populate("#locuszoom", data_sources, layout);

        $(".control-buttons .pan-left").on("click", function() {move(lzplot, -0.5);});
        $(".control-buttons .pan-right").on("click", function() {move(lzplot, 0.5);});
        $(".control-buttons .zoom-in").on("click", function() {zoom(lzplot, 1/1.5);});
        $(".control-buttons .zoom-out").on("click", function() {zoom(lzplot, 1.5);});
    });

});


