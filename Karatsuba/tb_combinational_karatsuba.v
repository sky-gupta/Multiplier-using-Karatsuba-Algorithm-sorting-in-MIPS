`timescale 1ns/1ps
`include "combinational_karatsuba.v"
module tb_combinational_karatsuba;

parameter N = 16;

// declare your signals as reg or wire
reg [N-1:0] X,Y;
wire c;
wire [2*N-1:0] Z;
wire [N-1:0] W;
initial begin

// write the stimuli conditions
assign X=16'b1111111111111111;
assign Y=16'b1111111111111111;

end

karatsuba_16 dut (.X(X), .Y(Y), .Z(Z));
// karatsuba_8 dut (.X(X), .Y(Y), .Z(Z));
// karatsuba_4 dut (.X(X), .Y(Y), .Z(Z));
// karatsuba_2 dut (.X(X), .Y(Y), .Z(Z));
// sub #(.N(4)) dut (.a(X),.b(Y),.cin(1'b1),.S(W),.cout(c));

initial begin
    // $dumpfile("combinational_karatsuba.vcd");
    // $dumpvars(0, tb_combinational_karatsuba);
    $monitor("X:%b,Y=%b,Z=%b",X,Y,Z);
end

endmodule
