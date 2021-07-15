classdef resonatorLC < handle
   
    properties
        fs
        oms
        C0
        Rm
        Lm
        Cm
        Qs
        kt2
        Rs
        Rp
        Nmodes
    end
    
    methods
        
        function obj = loadProperties(obj, fs, C0, Qs, kt2, Rs, Rp)
        
            obj.fs = fs;
            obj.C0 = C0;
            obj.Qs = Qs;
            obj.kt2 = kt2;
            obj.Rs = Rs;
            obj.Rp = Rp;
            obj.Nmodes = length(fs);
        
        end
        
        function obj = init(obj)
            
            obj.oms = 2*pi*obj.fs;
            obj.Rm = pi^2/8*1./(obj.oms.*obj.C0.*obj.kt2.*obj.Qs);
            obj.Lm = pi^2/8*1./(obj.oms.^2.*obj.C0.*obj.kt2);
            obj.Cm = 8/(pi^2)*obj.C0.*obj.kt2;
           
        end
        
        function Y = buildResponse(obj, f)
        
            om = 2*pi*f;
            Y0_i = 1i*om*obj.C0;
            Z0 = 1./Y0_i + obj.Rp;
            Y0 = 1./Z0;
            YBVD_all = zeros(1, length(f));
            
            for i = 1:obj.Nmodes
                
                YBVD_all = YBVD_all + 1./(obj.Rm(i) + 1i*om*obj.Lm(i) + 1./(1i*om*obj.Cm(i)));
                
            end
            
            YBVD_i = Y0 + YBVD_all;
            
            ZBVD = 1./YBVD_i + obj.Rs;
            Y = 1./ZBVD;
            
        end
        
    end
    
end