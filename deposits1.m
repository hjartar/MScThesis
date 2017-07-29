function y = deposits1(d)
global eta w0 a rd liq ws

%remember to replace "^" with ".^"
% fun = a*rd*((d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1))).^(1 - a)*(w0 + d*rd).^(a - 1) + (eta*(w0 + d*rd).^a*(d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1) - 1)*(a - 1)*(((w0 - d).^((eta - 1)/eta - 1)*(eta - 1))/eta - (d.^((eta - 1)/eta - 1)*liq*(eta - 1))/eta))/(((d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1))).^a*(eta - 1));
fun = a*((d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1))).^(1 - a)*(w0 + d*(rd - 1)).^(a - 1)*(rd - 1) + (eta*(w0 + d*(rd - 1)).^a*(d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1) - 1)*(a - 1)*(((w0 - d).^((eta - 1)/eta - 1)*(eta - 1))/eta - (d.^((eta - 1)/eta - 1)*liq*(eta - 1))/eta))/(((d.^((eta - 1)/eta)*liq + (w0 - d).^((eta - 1)/eta)).^(eta/(eta - 1))).^a*(eta - 1))

y=abs(fun);
