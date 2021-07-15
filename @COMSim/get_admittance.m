function [f,y]=get_admittance(obj)
    
    model = obj.model;
    f = mphglobal(model,'freq','dataset','dset7');
    y = mphglobal(model,'es2.I0_1/es2.V0_1','dataset','dset7');
    
end