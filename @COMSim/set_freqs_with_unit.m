function set_freqs_with_unit(obj, study, f, unit)

    model = obj.model;
 
    model.study(study).feature('freq').set('plist', f);
    model.study(study).feature('freq').set('punit', unit);

end