y1_temp = y1_explicit(end) + (y1_explicit(end) * exp(-x_explicit(end)^2) + x_explicit(end) * y2_explicit(end)) * h;
y2_temp = y2_explicit(end) + (3 * x_explicit(end) - y1_explicit(end) + 2 * y2_explicit(end)) * h;

y1_explicit_half = y1_explicit(end) + (y1_explicit(end) * exp(-x_explicit(end)^2) + x_explicit(end) * y2_explicit(end)) * (h / 2);
y2_explicit_half = y2_explicit(end) + (3 * x_explicit(end) - y1_explicit(end) + 2 * y2_explicit(end)) * (h / 2);

y1_explicit_double = y1_explicit_half + (y1_explicit_half * exp(-(x_explicit(end) + h / 2)^2) + (x_explicit(end) + h / 2) * y2_explicit_half) * (h / 2);
y2_explicit_double = y2_explicit_half + (3 * (x_explicit(end) + h / 2) - y1_explicit_half + 2 * y2_explicit_half) * (h / 2);