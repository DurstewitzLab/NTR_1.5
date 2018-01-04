function product=multinom_kernel(data_ref,order,data_val,beta,ganma)
%Private, auxiliar function.
%Computes the scalar product of any two vectors of a high-dimensional, feature space.
%High-dimensional space is defined by all possible products between
%original axes (i.e. inhomogenoues polynomial kernel).
%_________________________________________________________________________
%Inputs: 
%   beta, ganmma (1x1) = Constants of little relevance in general. For example, 
%                     expansion order O=2, if one wants to remove the multiplicity in the 
%                     products, choose beta=order.^(1/(order-1)); gannma=order.^(order/(1-order));
%
%   order (1x1) = Polynomial order of activity products.
%
%   data_ref (1 x n_dimensions (e.g. neurons)) = Must contain activity values over
%                                                     time.
%   data_val (1 x n_dimensions (e.g. neurons)) = Same for the validation data which 
%                                     will be then projected into the reference
%                                     eignevectors (please see comments in the 
%                                     code).
%Outputs: 
%   
%   product (1x1)= Value of the scalar product.
%_________________________________________________________________________
    
if nargin<4
    beta=1; ganma=1;
end
if nargin<3
    data_val=data;
end
product=ganma*(   ((1+ beta*data_val*data_ref')^order) - 1     );