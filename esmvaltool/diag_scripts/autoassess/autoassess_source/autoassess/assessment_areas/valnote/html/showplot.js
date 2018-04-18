function updateFields() {

fields = new Array(20);

var j = -1;

for (var i=0; i<nplot; i++) {
  j++;
  fields[j] = plots[i].field;
}

fields = fields.slice(0,j+1).unique();

var nfield = fields.length;


field = document.getElementById("field");

field.length = nfield;
for (var i=0; i<nfield; i++) {
    field.options[i].text = fields[i];
    field.options[i].value = fields[i];
}

}




function updateSeasons(element) {

var field = element.form.field.value;

seasons = new Array(5);
var j = -1;

for (var i=0; i<nplot; i++) {
    if (plots[i].field == field) {
	j++;
	seasons[j] = plots[i].season;
    }
}


seasons = seasons.slice(0,j+1).unique();
var nseas = seasons.length;

element.form.season.length = nseas;
for (var i=0; i<nseas; i++) {
    element.form.season.options[i].text = seasons[i];
    element.form.season.options[i].value = seasons[i];
}
if (nseas == 1) {
    element.form.season.disabled = true;
    /* update Obs menu - as won't be clicking season menu */
    updateObs(element);
} else {
    element.form.season.disabled = false;
}

}



function updateObs(element) {

var field = element.form.field.value;
var season = element.form.season.value;

obs = new Array(5);
var j = -1;

for (var i=0; i<nplot; i++) {

    if ((plots[i].field == field) &&
	(plots[i].season == season)) {
	j++;
	obs[j] = plots[i].obs;
    }
}


obs = obs.slice(0,j+1).unique();
var nobs = obs.length;

element.form.obs.length = nobs;
for (var i=0; i<nobs; i++) {
     element.form.obs.options[i].text = obs[i];
     element.form.obs.options[i].value = obs[i];
}

if (nobs == 1) {
    element.form.obs.disabled = true;
    /* update Proc menu - as won't be clicking obs menu */
    updateProc(element);
} else {
    element.form.obs.disabled = false;
}


}


function updateProc(element) {

var field = element.form.field.value;
var season = element.form.season.value;
var obs = element.form.obs.value;

proc = new Array(5);
var j = -1;

for (var i=0; i<nplot; i++) {

    if ((plots[i].field == field) &&
	(plots[i].season == season) &&
	(plots[i].obs == obs)) {
	j++;
	proc[j] = plots[i].proc;
    }
}


proc = proc.slice(0,j+1).unique();
var nproc = proc.length;

element.form.proc.length = nproc;
for (var i=0; i<nproc; i++) {
     element.form.proc.options[i].text = proc[i];
     element.form.proc.options[i].value = proc[i];
}

if (nproc == 1) {
    element.form.proc.disabled = true;
} else {
    element.form.proc.disabled = false;
}


}



function showplot(element) {

var fname = filename(element);
var html = "<img src='"+fname+"'/>";

/* todo: record previous filenames, and add back/forward buttons */

document.getElementById("plot_box").innerHTML=html;
}




function filename(element) {

var field = element.form.field.value;
var season = element.form.season.value;
var obs = element.form.obs.value;
var proc = element.form.proc.value;



for (var i=0; i<nplot; i++) {
    if ((plots[i].field == field) &&
	(plots[i].season == season) &&
	(plots[i].obs == obs) &&
	(plots[i].proc == proc)) {
	var file = plots[i].file;
	break;
    }
}

return file;

}


function Plot(file, field, season, obs, proc) {
this.file = file;
this.field = field;
this.season = season;
this.obs = obs;
this.proc = proc;
}


// New method for arrays - return array of unique elements
Array.prototype.unique =
function() {
  var u = [];
  var ulen = 0;
  var l = this.length;
  for(var i=0; i<l; i++) {
    first = true;
    for(var j=0; j<ulen; j++) {
      // is the first instance of this element?
      // i.e. have we put in in 'a' already?
      if (this[i] == u[j]) {
        first = false;
        break;
      }
    }
    if (first) {
      ulen = u.push(this[i]);
    }
  }
  return u;
};
