//UI Bindings using the Data_manager type and D3.js.

Plot = function(element, data, display, options, viewscope){
	this.element = element;
	this.data = data;
	this.display = new display(element, data, viewscope, options);
	this.display.draw();
	this.redraw = function(){
		if (typeof(this.display.redraw) == "function"){
		this.display.redraw();
	}
	}
}

UI_Manager = function(viewscope){

	this.viewscope = viewscope;
	this.viewscope.subscribe(this);

	this.plots = [];

	//Special case binding
	this.bind_view_scope = function(element){
		this.plots.push(new Plot(element, this.viewscope, ViewSlider, {}));
	}

	this.notify = function(type, data){
		for (var pid in this.plots){
			this.plots[pid].redraw();
		}
	}
	
	this.bind = function(element, data, display, options){
		options = (options===undefined) ? {} : options;
		this.plots.push(new Plot(element, data, display, options,this.viewscope));
	}
}

D3AnnotationPlot = function(element, data, viewscope, options){
	//Displays a list of annotations on a line.
	this.element = element;
	this.data = data;
	this.viewscope = viewscope;
	this.options = options;
	var self = this;
	this.draw = function(){
		PH = this.data.annotations.length * 5 + 40;
		self.x_scale = d3.scale.linear().domain([self.viewscope.min, self.viewscope.max]).range([0, PLOT_WIDTH - PLOT_MARGIN * 2]);

		self.chart = d3.select(element).append("svg:svg")
            .attr("viewBox", "0 0 " + PLOT_WIDTH + " " + PH)
            .attr("preserveAspectRatio", "xMidYMid meet")
            .attr("class", "figure cscheme_white");
        var x_axis = d3.svg.axis().scale(self.x_scale).ticks(10).orient("bottom");

        self.chart.append("svg:g")
            .attr("class", "axis y_axis cscheme cscheme_white")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + (PH - PLOT_MARGIN) + ")")
            .call(x_axis);
        for (var aid in this.data.annotations){
            a = this.data.annotations[aid];
            c = "white";
            st = JSON.stringify(a.summary);
            if (st.search("phyre2") != -1){
            	c = "NavajoWhite";
            }
            if (st.search("Coil") != -1){
            	c = "PaleGreen";
            }
            if (st.search("Region of a membrane") != -1){
            	c = "PaleTurquoise";
            }
            if (st.search("SSF") != -1){
            	c = "Pink";
            }
            if (st.search("Pfam") != -1){
            	c = "Plum";
            }
            self.chart.append("svg:rect")
            .attr("x", PLOT_MARGIN + self.x_scale(a.start))
            .attr("y", aid * 5 )
            .attr("width", self.x_scale(Math.min(1000,a.end) - a.start))
            .attr("height", 5)
            .attr("title",a.title)
            .attr("fill", c);
        }
	}
}

AnnotationList = function(element, data, viewscope, options){
	//Displays a simple annotation list on screen. Data should be an annotation list.
	this.container = $(element);
	this.data = data;
	this.options = options;
	this.viewscope = viewscope;
	var self = this;					

	this.draw = function(){
		this.container.append(nel("h2", "All Domains"));
		for (var source in this.data.sources){
			this.container.append(nel("h3", source))
			for (var annotation in this.data.sources[source]){
				new_element = annotation_as_html(this.data.annotations[this.data.sources[source][annotation]])
				this.container.append(new_element);
				new_element.find("div").css("max-height", new_element.find("div").height() + "px")
				new_element.addClass("collapsed");
				new_element.find(".btn").click(function(e){
				});
			}
		}
	}
}

WriteMeta = function(element, data, viewscope, options){
	this.container = $(element);
	this.data = data;
	this.options = options;

	this.draw = function(){
		this.container.text(data[this.options.key]);
	}
}

PLOT_WIDTH = 500;
PLOT_HEIGHT = 120;
PLOT_MARGIN = 25;
DEFAULT_CS = ["cscheme_lightyellow", "cscheme_lightgreen", "cscheme_lightblue"];

D3SeriesGroupOverview = function(element, data, viewscope, options) {

    this.container = element;
    this.data = data;
    this.options = options;
    this.viewscope = viewscope;

    var self = this;

    this.redraw = function() {
        self.data.start = this.viewscope.start;
        self.data.end = this.viewscope.end;
        s = [{
            position: self.data.start,
            value: self.ymax
        }, {
            position: self.data.end,
            value: self.ymax
        }];
        select_area = self.chart.select(".select_area");
        select_area.datum(s);
        select_area.transition(0)
            .attr("d", self.area_plot);
    };
    this.draw = function() {
        //Generate values to plot
        self.series_names = Object.keys(self.data.series);
        self.ymin = Infinity;
        self.ymax = -Infinity;
        self.min = self.data.series[self.series_names[0]].series.min;
        self.max = self.data.series[self.series_names[0]].series.max;
        self.start = self.data.series[self.series_names[0]].series.start;
        self.end = self.data.series[self.series_names[0]].series.end;
        var allvalues = [];
        for (var series_id in self.series_names){
        	this_series = self.data.series[self.series_names[series_id]];
        	self.ytlower = this_series.ymin;
        	self.ytupper = this_series.ymax;
        	self.ymin = Math.min(self.ymin, _.min(this_series.series.d()));
        	self.ymax = Math.max(self.ymax, _.max(this_series.series.d()));
        	allvalues[self.series_names[series_id]] = transpose({
            	position: this_series.series.p(),
            	value: this_series.series.d(),
            	value_smooth: gauss_smooth(this_series.series.d(), this_series.gauss)
        	});

        }

        if (options.round) {
            self.ymin = Math.floor(self.ymin);
            self.ymax = Math.ceil(self.ymax);
        }

        //Set up geometry generators
        self.area_plot = d3.svg.area()
            .interpolate("step")
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y1(function(d) {
                return self.y_scale(d.value);
            })
            .y0(function() {
                return self.y_scale(self.ymin);
            });

        var th_area_plot = d3.svg.area()
            .interpolate("step")
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y1(function(d) {
                if (d.value_smooth > self.ytlower && d.value_smooth < self.ytupper) {
                    return self.y_scale(d.value);
                } else {
                    return self.y_scale(self.ymin);
                }
            })
            .y0(function() {
                return self.y_scale(self.ymin);
            });

        self.line_plot = d3.svg.line()
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y(function(d) {
                return self.y_scale(d.value);
            });

        //Set up scales
        self.x_scale = d3.scale.linear().domain([self.min, self.max]).range([0, PLOT_WIDTH - PLOT_MARGIN * 2]);
        self.y_scale = d3.scale.linear().domain([self.ymin, self.ymax]).range([PLOT_HEIGHT - PLOT_MARGIN * 2, 0]);

        //Create chart
        self.chart = d3.select(element).append("svg:svg")
            .attr("viewBox", "0 0 " + PLOT_WIDTH + " " + PLOT_HEIGHT)
            .attr("preserveAspectRatio", "xMidYMid meet")
            .attr("class", "figure");
            csid = -1;
         for (var v in allvalues){
         	csid += 1;
         	values = allvalues[v];
        //Add area and line representations
        var legend = self.chart.append("svg:g")
             .attr("class", "legend " + self.options.colours[csid]);
        legend.append("svg:text")
             .attr("transform", "translate(" + (PLOT_MARGIN*2 + (50*csid) ) + "," + (PLOT_MARGIN) + ")")
             .attr("x2", 10)
             .attr("class", "legend cf1 " + self.options.colours[csid])
             .attr("text-anchor", "middle")
             .text(v);

         self.chart.append("svg:g")
             .attr("style", "mask: url(#selection_mask)")
             .append("svg:path")
             .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
             .datum(values)
            .attr("d", th_area_plot)
             .attr("class", "fig_area cf1 " + self.options.colours[csid]);

        self.chart.append("svg:g")
             .attr("style", "mask: url(#selection_mask)")
            .append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
           .datum(values)
           .attr("d", self.line_plot)
            .attr("class", "fig_line cs1 " + self.options.colours[csid]);
    	};

        //Add selection mask
        var selection_mask = self.chart.append("svg:defs")
            .append("mask")
            .attr("id", "selection_mask");
        selection_mask.append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum([{
                position: self.min,
                value: this.ymax
            }, {
                position: self.max,
                value: this.ymax
            }])
            .attr("d", self.area_plot)
            .attr("class", "select_area_negative");
        selection_mask.append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum([{
                position: self.start,
                value: this.ymax
            }, {
                position: self.end,
                value: this.ymax
            }])
            .attr("d", self.area_plot)
            .attr("class", "select_area");

        //Plot upper and lower thresholds
        var th = self.chart.append("svg:g")
            .attr("class", "threshold cs1  cscheme_" + self.options.cscheme)
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")");
        th.append("svg:path")
            .datum([{
                position: self.start,
                value: self.ytlower
            }, {
                position: self.end,
                value: self.ytlower
            }])
            .attr("d", self.line_plot);
        th.append("svg:path")
            .datum([{
                position: self.start,
                value: self.ytupper
            }, {
                position: self.end,
                value: self.ytupper
            }])
            .attr("d", self.line_plot);

        //Plot Axes
        var x_axis = d3.svg.axis().scale(self.x_scale).ticks(10).orient("bottom");

        self.chart.append("svg:g")
            .attr("class", "axis y_axis cscheme_" + self.options.cscheme)
            .attr("transform", "translate(" + PLOT_MARGIN + "," + (PLOT_HEIGHT - PLOT_MARGIN) + ")")
            .call(x_axis);

        self.y_axis_tick_values = [self.ymin, self.ymax, self.ytupper, self.ytlower];
        var y_axis = d3.svg.axis().scale(self.y_scale).tickValues(self.y_axis_tick_values).orient("left");

        self.chart.append("svg:g")
             .attr("class", "axis x_axis cscheme_" + self.options.cscheme)
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .call(y_axis);
        this.redraw();
    };
};


D3SingleSeriesOverview = function(element, data, viewscope, options) {

    this.container = element;
    this.data = data;
    this.options = options;
    this.viewscope = viewscope;

    var self = this;

    this.redraw = function() {
        self.data.series.start = this.viewscope.start;
        self.data.series.end = this.viewscope.end;
        s = [{
            position: self.data.series.start,
            value: self.ymax
        }, {
            position: self.data.series.end,
            value: self.ymax
        }];
        select_area = self.chart.select(".select_area");
        select_area.datum(s);
        select_area.transition()
            .attr("d", self.area_plot);
    };
    this.draw = function() {
        //Set up initial values
        self.ymin = _.min(this.data.series.d());
        self.ymax = _.max(this.data.series.d());
        self.ytlower = this.data.ymin;
        self.ytupper = this.data.ymax;

        if (options.round) {
            self.ymin = Math.floor(self.ymin);
            self.ymax = Math.ceil(self.ymax);
        }

        //Set up geometry generators
        self.area_plot = d3.svg.area()
            .interpolate("step")
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y1(function(d) {
                return self.y_scale(d.value);
            })
            .y0(function() {
                return self.y_scale(self.ymin);
            });

        var th_area_plot = d3.svg.area()
            .interpolate("step")
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y1(function(d) {
                if (d.value_smooth > self.ytlower && d.value_smooth < self.ytupper) {
                    return self.y_scale(d.value);
                } else {
                    return self.y_scale(self.ymin);
                }
            })
            .y0(function() {
                return self.y_scale(self.ymin);
            });

        self.line_plot = d3.svg.line()
            .x(function(d) {
                return self.x_scale(d.position);
            })
            .y(function(d) {
                return self.y_scale(d.value);
            });

        //Set up scales
        self.x_scale = d3.scale.linear().domain([self.data.series.min, self.data.series.max]).range([0, PLOT_WIDTH - PLOT_MARGIN * 2]);
        self.y_scale = d3.scale.linear().domain([self.ymin, self.ymax]).range([PLOT_HEIGHT - PLOT_MARGIN * 2, 0]);

        //Generate values to plot
        var values = transpose({
            position: this.data.series.p(),
            value: this.data.series.d(),
            value_smooth: gauss_smooth(this.data.series.d(), this.data.gauss)
        });

        //Create chart
        self.chart = d3.select(element).append("svg:svg")
            .attr("viewBox", "0 0 " + PLOT_WIDTH + " " + PLOT_HEIGHT)
            .attr("preserveAspectRatio", "xMidYMid meet")
            .attr("class", "figure cscheme_" + self.options.cscheme);

        //Add area and line representations
        self.chart.append("svg:g")
            .attr("style", "mask: url(#selection_mask)")
            .append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum(values)
            .attr("d", th_area_plot)
            .attr("class", "fig_area cf1");

        self.chart.append("svg:g")
            .attr("style", "mask: url(#selection_mask)")
            .append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum(values)
            .attr("d", self.line_plot)
            .attr("class", "fig_line cs1");


        //Add selection mask
        var selection_mask = self.chart.append("svg:defs")
            .append("mask")
            .attr("id", "selection_mask");
        selection_mask.append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum([{
                position: this.data.series.min,
                value: this.ymax
            }, {
                position: this.data.series.max,
                value: this.ymax
            }])
            .attr("d", self.area_plot)
            .attr("class", "select_area_negative");
        selection_mask.append("svg:path")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .datum([{
                position: this.data.series.start,
                value: this.ymax
            }, {
                position: this.data.series.end,
                value: this.ymax
            }])
            .attr("d", self.area_plot)
            .attr("class", "select_area");

        //Plot upper and lower thresholds
        var th = self.chart.append("svg:g")
            .attr("class", "threshold cs1")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")");
        th.append("svg:path")
            .datum([{
                position: this.data.series.start,
                value: self.ytlower
            }, {
                position: this.data.series.end,
                value: self.ytlower
            }])
            .attr("d", self.line_plot);
        th.append("svg:path")
            .datum([{
                position: this.data.series.start,
                value: self.ytupper
            }, {
                position: this.data.series.end,
                value: self.ytupper
            }])
            .attr("d", self.line_plot);

        //Plot Axes
        var x_axis = d3.svg.axis().scale(self.x_scale).ticks(10).orient("bottom");

        self.chart.append("svg:g")
            .attr("class", "axis y_axis")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + (PLOT_HEIGHT - PLOT_MARGIN) + ")")
            .call(x_axis);

        self.y_axis_tick_values = [self.ymin, self.ymax, self.ytupper, self.ytlower];
        var y_axis = d3.svg.axis().scale(self.y_scale).tickValues(self.y_axis_tick_values).orient("left");

        self.chart.append("svg:g")
            .attr("class", "axis x_axis")
            .attr("transform", "translate(" + PLOT_MARGIN + "," + PLOT_MARGIN + ")")
            .call(y_axis);
        this.redraw();
    };
};

PLOT_PER_LINE = 30;

D3SecondaryStructure = function(element, data, viewscope, options) {

    this.container = element;
    this.data = data;
    this.options = options;
    this.viewscope = viewscope;

    var self = this;


    this.draw = function(){
        this.root = d3.select(this.container);
        //Number of bases to plot:
        n = (this.viewscope.start - this.viewscope.min) + 1;
        //Each base gets a box of width PLOT_WIDTH/PLOT_PER_LINE
        // and height = width * 2
        //N bases will need ceil(N/PLOT_PER_LINE) lines
        // so height = PLOT_WIDTH*2/PLOT_PER_LINE * ceil(N/PLOT_PER_LINE)
        plot_height = PLOT_WIDTH*2/PLOT_PER_LINE * Math.ceil(n/PLOT_PER_LINE);

        self.chart = this.root.append("svg:svg")
            .attr("viewbox", "0 0 " + PLOT_WIDTH + " 0")
            .attr("preserveAspectRatio", "xMidYMid meet")
            .attr("class", "figure figure_secondary-structure");
        for (var pos = this.viewscope.min; pos < this.viewscope.max; pos++){
            var aa = this.data.sequence.series.data[pos];
            var characters = []
            characters.push(["c",this.data.secondary_structure.series.Coil.series.data[pos]]);
            characters.push(["s",this.data.secondary_structure.series.Strand.series.data[pos]]);
            characters.push(["h",this.data.secondary_structure.series.Helix.series.data[pos]]);
            characters.sort(function(a,b){return b[1]-a[1]});
            var base_plot = self.chart.append("g")
                .attr("class", "ss_base base_"+pos+" IUPAC_"+aa)
            base_plot.append("text").text(aa)
                .attr("class", "aa_code").attr("y",-10).attr("text-anchor", "center");;
            postext = "";
            if (pos % 5 == 0){
                postext = "|";
            }
            if (pos % 10 == 0){
                postext = pos;
            }
            base_plot.append("text").text(postext).attr("class", "aa_pos").attr("text-anchor", "center");
            if (characters[0][0] == "c"){
                base_plot.append("rect")
                    .attr("y",-35)
                    .attr("width", 1+PLOT_WIDTH/PLOT_PER_LINE)
                    .attr("height", 2)
                    .style("fill", "black");
            }
            if (characters[0][0] == "s"){
                base_plot.append("rect")
                    .attr("y",-39)
                    .attr("width", 1+PLOT_WIDTH/PLOT_PER_LINE)
                    .attr("height", 10)
                    .style("fill", "blue");
            }
            if (characters[0][0] == "h"){
                base_plot.append("rect")
                    .attr("y",-39)
                    .attr("width", 1+PLOT_WIDTH/PLOT_PER_LINE)
                    .attr("height", 10)
                    .style("fill", "green");
            }
        }
        self.redraw();
    }
    this.redraw = function(){
        //Number of bases to plot:
        n = (this.viewscope.end - this.viewscope.start) + 1;
        //Each base gets a box of width PLOT_WIDTH/PLOT_PER_LINE
        // and height = width * 2
        //N bases will need ceil(N/PLOT_PER_LINE) lines
        // so height = PLOT_WIDTH*2/PLOT_PER_LINE * ceil(N/PLOT_PER_LINE)
        plot_height = PLOT_WIDTH*3/PLOT_PER_LINE * Math.ceil(n/PLOT_PER_LINE);

        self.chart.attr("viewbox", "0 0 " + PLOT_WIDTH + " " + plot_height);
        self.chart.attr("height", plot_height)

        self.data.sequence.series.start = self.viewscope.start;
        self.data.sequence.series.end = self.viewscope.end;
        self.data.secondary_structure.start = self.viewscope.start;
        self.data.secondary_structure.end = self.viewscope.end;
        for (var pos = self.viewscope.min; pos < self.viewscope.max; pos++){
            offsetpos = pos - self.viewscope.start;
            var xpos = (PLOT_WIDTH/PLOT_PER_LINE) * (offsetpos % PLOT_PER_LINE);
            var ypos = (PLOT_WIDTH*3/PLOT_PER_LINE) * (1 + Math.floor(offsetpos / PLOT_PER_LINE));
            var base = self.chart.select(".base_"+pos);
            if (pos < self.viewscope.start || pos > self.viewscope.end){
                base.style("opacity",0);
            } else {
                base.style("opacity",1);
                base.attr("transform","translate("+xpos+","+ypos+")");
            }
        }
    }

}

MatchVisibility =  function(element, data, viewscope, options){
	this.container = $(element);
	this.data = data;
	this.options = options;

	this.draw = function(){
		display = "none";
		if (this.data[this.options["key"]]){
			display = "";
		}
		this.container.css("display", display);
	}
}

ViewSlider = function(element, data, viewscope, options){
	this.container = $(element);
	this.data = data;
	this.options = options;
	var self = this;
	this.draw = function(){
		this.container.attr("data-start",this.data.start)
		.attr("data-end",this.data.end)
		.attr("data-min",this.data.min)
		.attr("data-max",this.data.max)
		.attr("data-int", true);

		this.ui_element = slider(this.container);
		this.ui_element.on("update", function(e){
			self.data.start = self.ui_element.data("start");
			self.data.end = self.ui_element.data("end");
			self.data.update();
		});
		}
	}

annotation_as_html = function(annotation){
	var html = nel("section");
	var title = nel("h4", "<span>&#x276F;</span><a class = 'btn btn-xs btn-danger' title='"+annotation.title+"' href='#"+annotation.title + "-" + annotation.hash+"'>(" + annotation.start + " - " + annotation.end+")</a> " + annotation.title, "annotation_title");
	var properties_list = nel("dl", "", "annotation_properties");
	properties_list.append("<a class = 'btn btn-sm pull-right btn-danger' title='"+annotation.title+"' href='#"+annotation.title + "-" + annotation.hash+"'>View in detail</a>")
	properties_list.append(nel("dt", "Range"));
	properties_list.append(nel("dd", annotation.start + " - " + annotation.end));
	properties_list.append(nel("dt", "Source"));
	properties_list.append(nel("dd", annotation.source));
	for (property in annotation.summary){
		properties_list.append(nel("dt", property));
		properties_list.append(nel("dd", annotation.summary[property]));
	};
	title.click(function(e){
			sect = $(e.target).parents("section");
				sect.siblings("section").addClass("collapsed");
				sect.toggleClass("collapsed");
			});
	html.append(title).append(nel("div").append(properties_list));
	html.find(".btn").attr("data-start", annotation.start);
	html.find(".btn").attr("data-end", annotation.end);
	return html;
}

nel = function(eltag, elcontent, elclass){
	elclass = (elclass === undefined) ? "" : elclass;
	elcontent = (elcontent === undefined) ? "" : elcontent;
	return $('<'+eltag+' class="'+elclass+'">'+elcontent+'</'+eltag+">");
}

makeColourScheme = function(name, colour){
cs_rules = "<style>.cscheme_NAME.cs1,.cscheme_NAME .domain, .cscheme_NAME .cs1, .cscheme_NAME.axis line, .cscheme_NAME.axis path .cscheme_NAME .axis line, .cscheme_NAME .axis path{stroke: COLOUR;} .cscheme_NAME.cf1, .cscheme_NAME .cf1, .cscheme_NAME text, .cscheme_NAME .fig_area{fill: COLOUR;}</style>";
cs_rules = cs_rules.replace(/COLOUR/g, colour);
cs_rules = cs_rules.replace(/NAME/g, name);

$('html > head').append($(cs_rules));
}