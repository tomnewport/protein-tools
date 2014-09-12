//Data manager.
//Requeres: underscore.js
//          jquery
//          numeric.js

function tic(){
	tictime = new Date();
}

function toc(){
	var toctime = new Date();
	return ((toctime - tictime)/1000);
}

function fill(width, value){
	var d = [];
	for (var x = 0; x<width; x++){
		d.push(value);
	}
	return d;
}

function gauss(x,y,sigmaSqr){
	//from http://www.embege.com/gauss/
	var pow =  Math.pow(Math.E, -( (x*x+y*y) / (2*sigmaSqr) ));
	var g = ( 1/(Math.sqrt(2*Math.PI*sigmaSqr)) ) * pow; 
	return g;
}

function kernel(type, size){
	var adds_to_one = true;
  var k = [];
	//Generates a 1D kernel for convolution
	if (type == "uniform" || type === undefined){
		//Uniform density
    k = fill(size, 1);
	}
	if (type == "gauss"){
		//Gaussian falloff
		var cutoff = 0.0000001;
		k = [gauss(0,0,size)];
		var x = 0;
		while (k[0] > cutoff){
			x++;
			var g = gauss(x,0,size);
			k = Array.concat(g,k,g);
		}
	}
	if (type == "sobel"){
		//Sobel kernel
		adds_to_one = false;
		k = [0];
		for (var s = 1; s<=size; s++){
			k = Array.concat(s,k,-s);
		}
	}
	//Normalise kernel if it adds to one

	k = (adds_to_one)? numeric["/"](k, numeric.sum(k)) : k;
	return k;
}

function apply_boundary(array, method, size){
	if (method == "uniform"){
		var left = fill(Math.ceil(size/2), array[0]);
		var right = fill(Math.floor(size/2), array[array.length-1]);
		var array = left.concat(array).concat(right);
	}
	return array;
}

function convolve(kernel, array, boundary){
	// 1D kernel convolution with reasonable performance.
	var result = [];
	array = apply_boundary(array, boundary, kernel.length - 1);
	for (var x = 0; x<=array.length - kernel.length; x++){
		var s = array.slice(x,x+kernel.length);
		result[x] = numeric.sum( numeric["*"](s, kernel));
	}
	return result;
}

function gauss_smooth(data, size){
	return convolve(kernel("gauss",size), data, "uniform");
}

function in_range(data, min, max){
	var above = numeric["<"](data, max);
	var below = numeric[">"](data, min);
	return {"above":above, "below":below, "within": numeric["&&"](above, below)};
}

function compress(array, offset){
	//Function accepts an array. No element in the array should be undefined.
	//Returned is an object of the following scheme:
	//	root
	//		state
	//			start : end
	//			start : end
	//		state
	//			start : end

	array.push(undefined);
	var pos = (offset === undefined)? 0 : offset;
	var state = array[0];
	var labels = {};
	var start = pos;
	for (var element in array){
		if (array[element] != state){
			if (!(state in labels)){
				labels[state] = {};
			}
			labels[state][start] = pos - 1;

			state = array[element];
			start = pos;
		}
		pos += 1;
	}
	return labels;
}

function Series(array, options){
	//A cached slicable array.
	this.data = array;
	options = (options===undefined)? {} : options;
	options.min = (options.min===undefined)? 0 : eval(options.min);
	options.max = (options.max===undefined)? array.length : eval(options.max);
	options.start = (options.start===undefined)? options.min : eval(options.start);
	options.end = (options.end===undefined)? options.max : eval(options.end);
	//start and end are INCLUSIVE
	this.start = options.start;
	this.end = options.end;
	this.max = options.max;
	this.min = options.min;

	//Check for changes to the object and update cache if needed. 
	//Using cache makes operations several orders of magnitude faster.
	//update_cache(true) should be called every time this.data, this.min and this.max are changed.
	this.update_cache = function(force){
		if (force === true || this.old_start != this.start || this.old_end != this.end){
			//Cache is out of date or update has been forced
			this.old_start = this.start;
			this.old_end = this.end;
			this.cache_d = this.data.slice(this.start-this.min, this.end + 1);
			this.cache_p = _.range(this.start, this.end + 1);
			this.cache_o = {};
			for (var x in _.range(0,this.cache_d.length)){
				this.cache_o[ this.cache_p[x] ] = this.cache_d[x];
			}
		}
		return this;
	};

	//Sets the range of data
	this.set_range = function(start, end){
		this.start = start;
		this.end = end;
		return this;
	};
	//Returns sliced data
	//If force is true, updates the cache, otherwise check if start and end have changed.
	this.d = function(force){
		return this.update_cache(force).cache_d;
	};
	//Returns positions corresponding to data
	this.p = function(force){
		return this.update_cache(force).cache_p;
	};
	//Returns an object of position:data
	this.o = function(force){
		return this.update_cache(force).cache_o;
	};
	this.annotate = function(min, max, gauss, title){
		var opts = {"source":"generated", "title": title, "hash": "", "summary": {"Range":min + " - " + max, "Smoothing": gauss}};
		var ranges = in_range(gauss_smooth(this.d(), gauss), min, max);
		var labels = compress(ranges.within);
		var result = [];
		for (var start in labels.true){
			var newLabel = {};
			newLabel.start = start;
			newLabel.end = labels.true[start];
			newLabel = _.extend(newLabel, opts);
			result.push(new Annotation(newLabel));
		}
		return new Annotation_list(result);
	};

	this.update_cache(true);
	return this;
}

function Annotation(obj){
	// This is a single annotation of a genomic region
	if (_.isString(obj.summary)){
		obj.summary = JSON.parse(obj.summary);
	}
	this.start = obj.start;
	this.end = obj.end;
	this.title = obj.title;
	this.hash = obj.hash;
	this.source = obj.source;
	this.summary = obj.summary;
	return this;
}

function Annotation_list(list){
	// This is a cached and indexed list of annotations
	this.annotations = (list === undefined)? [] : list;
	this.tag_heirachy = [];
	var d = this.annotations;

	this.in_range = function(start, end){
		var result = [];
		for (var ind in this.annotations){
			if (this.annotations[ind].start > start || this.annotations[ind].end < end){
				result.push(this.annotations[ind]);
			}
		}
	};
	this.max = function(){
		return _.max(this.transposed.end);
	};
	this.min = function(){
		return _.min(this.transposed.start);
	};
	this.non_overlap = function(options){
		//TODO
	};
	this.by_source = function(source){
		var result = [];
		source_ids = this.sources[source];
		for (var id in source_ids){
			result.push(this.annotations[source_ids[id]]);
		}
		return result;
	}
	this.refresh = function(){
		//Allow manipulation of summary elements?
		this.sources = {};
		this.transposed = {};
		this.transposed.start = [];
		this.transposed.end = [];
		this.transposed.title = [];
		this.transposed.hash = [];
		this.transposed.summary = [];
		this.transposed.source = [];
		for (var ind in this.annotations){
			var annotation = this.annotations[ind];
			this.transposed.source[ind] = annotation.source;
			this.transposed.hash[ind] = annotation.hash;
			this.transposed.start[ind] = annotation.start;
			this.transposed.end[ind] = annotation.end;
			this.transposed.title[ind] = annotation.title;
			this.transposed.summary[ind] = annotation.summary;
			if (!(annotation.source in this.sources)){
				this.sources[annotation.source] = [];
			}
			this.sources[annotation.source].push(ind);
		}
	};

	this.refresh();
	return this;
}

function Annotated_series(series, ymin, ymax, gauss, title){
	this.ymin = ymin;
	this.ymax = ymax;
	this.gauss = gauss;
	this.series = series;
	this.annotations = this.series.annotate(ymin, ymax, gauss, title);
}

function Series_group(){
	this.series = {};
	var self = this;
	this.__defineSetter__("start", function(val){
        for (seriesname in self.series){
        	if (_.isObject(self.series[seriesname])){
        		self.series[seriesname].series.start = val;
        	}
        }
        this._start = val;
        return val;
    });
    this.__defineGetter__("start", function(val){
        return this._start;
    });

	this.__defineSetter__("end", function(val){
        for (seriesname in self.series){
        	if (_.isObject(self.series[seriesname])){
        		self.series[seriesname].series.end = val;
        	}
        }
        this._end = val;
        return val;
    });
    this.__defineGetter__("end", function(val){
        return this._end;
    });

	this.sum = false;
	this.parameterspace = false;
}

function transpose(in_unknown) {
    if (_.isArray(in_unknown)) {
        var in_array = in_unknown;
        var out_object = {};
        //Takes an array of objects and returns an object of arrays
        var keys = Object.keys(in_array[0]);
        for (var key_id in keys) {
            out_object[keys[key_id]] = [];
        }
        for (var item_id in in_array) {
            for (key_id in keys) {
                out_object[keys[key_id]][item_id] = in_array[item_id][keys[key_id]];
            }
        }
        return out_object;

    } else {
        var in_object = in_unknown;
        var out_array = [];
        //Takes an object of arrays and returns an array of objects
        var keys = Object.keys(in_object);
        for (var item_id in in_object[keys[0]]) {
            out_array[item_id] = {};
            for (var key_id in keys) {
                out_array[item_id][keys[key_id]] = in_object[keys[key_id]][item_id];
            }
        }
        return out_array;

    }
}

function Data_manager(){
	var self = this;
	this.json_schema = {};
	this.metadata = {};
	this.csv_files = {};
	this.min = 0;
	this.max = 0;
	this.series = {};
	this.series_groups = {};
	this.annotations = {};
	this.user_annotations = {};
	this.new_user_annotation = {};
	this.proceed_if_loaded = function(if_loaded){
		//Performs the function if_loaded
		loaded = true;
		for (var file_id in this.csv_files){
			if ( this.csv_files[file_id].loading === true ){
				loaded = false;
			}
		}
		if (this.metadata.loading === true){
			loaded = false;
		}
		if (loaded){
			if_loaded();
		}
	}
	this.from_json = function(filename){
		$.getJSON(filename, function(data){
				//Load all csv files
				self.json_schema = data.protein;
				for (var file_id in data.protein.csv_files){
					self.load_csv(data.protein.csv_files[file_id]);
				}
				//Load metadata
				self.load_meta(data.protein.metadata);
			});
	};
	this.get_series = function(identifier){
		idparts = identifier.split(".csv$");
		filename = idparts[0] + ".csv";
		colname = idparts[1];
		return this.csv_files[filename].as_object[colname];
	};
	this.load_meta = function(filename){
		//Loads and overwrites metadata
		this.metadata.loading = true;
		$.getJSON(filename, function(newmetadata){
			this.metadata = _.extend(self.metadata, newmetadata);
			self.metadata.loading = false;
			self.proceed_if_loaded(self.setup);
		});
	};
	this.load_csv = function(filename){
		self.csv_files[filename] = {};
		self.csv_files[filename].as_array = [];
		self.csv_files[filename].as_object = {};
		self.csv_files[filename].loading = true;
		d3.csv(filename, function(data){
			if (_.isArray(data)){
				self.csv_files[filename].as_array = data;
				self.csv_files[filename].as_object = transpose(data);
			} else {
				console.log("Warning - csv may not have been correctly loaded!");
				self.csv_files[filename].as_array.push(data);
			}
			self.csv_files[filename].loading = false;
			//If all files are loaded, set everything up:
			self.proceed_if_loaded(self.setup);
		});
	};
	this.setup = function(){
		//Function should be called once all files are loaded.

		//Set the start and end points for series.
		//ASSUMPTION: positions are a continuous series of consecutive ascending integers from start to finish
		series_pos = self.get_series(self.json_schema.position);
		self.min = _.min(series_pos);
		self.max = _.max(series_pos);
		//Set up series
		for (var series_id in self.json_schema.series){
			var series_schema = self.json_schema.series[series_id];
			var new_array = self.get_series(series_schema.file);
			var new_series = new Series(new_array, self);
			self.series[series_id] = new Annotated_series(new_series,series_schema.ymin,series_schema.ymax,series_schema.gauss,series_schema.annotation);
		}
		//Set up series groups
		for (var sgroup_id in self.json_schema.series_groups){
			var sgroup_schema = self.json_schema.series_groups[sgroup_id];
			var new_sgroup = new Series_group();
			for (var s_id in sgroup_schema.series){
				var new_array = self.get_series(sgroup_schema.series[s_id]);
				new_series = new Series(new_array, self);
				new_sgroup.series[s_id] = new Annotated_series(new_series, sgroup_schema.ymin, sgroup_schema.ymax, sgroup_schema.gauss, s_id);
			}
			new_sgroup.sum = sgroup_schema.sum;
			new_sgroup.parameterspace = sgroup_schema.parameterspace;
			self.series_groups[sgroup_id] = new_sgroup;
		}
		//Set up annotations
		var annotation_list = [];
		for (var annotation_id in self.json_schema.annotations){
			var new_annotations = self.csv_files[self.json_schema.annotations[annotation_id]].as_array;
			for (a_id in new_annotations){
				new_annotations[a_id] = new Annotation(new_annotations[a_id]);
			}
			annotation_list = annotation_list.concat(new_annotations);
		}
		a = annotation_list;
		self.annotations = new Annotation_list(annotation_list);
	self.ready();
	}
}

Viewscope = function(protein_data){
	this.protein_data = protein_data;
	this.start = protein_data.min;
	this.end = protein_data.max;
	this.min = protein_data.min;
	this.max = protein_data.max;
	this.view_domain = "";
	this.view_preposition = "";
	this.detail_view = false;
	this.subscribers = []
	this.subscribe = function(newSubscriber){
		this.subscribers.push(newSubscriber);
	}
	this.update = function(){
		for (var subscriber in this.subscribers){
			this.subscribers[subscriber].notify("viewscope",this);
		}
	}
}