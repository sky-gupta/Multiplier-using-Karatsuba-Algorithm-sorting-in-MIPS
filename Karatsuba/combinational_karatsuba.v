//ADDER
module half_adder(a, b, S, cout);
input a,b;
output S,cout;

assign S=a^b;
assign cout = a&b;
endmodule


module full_adder(a, b, cin, S, cout);
input a,b,cin;
output S,cout;
wire s1,c1,c2;

half_adder h1(a,b,s1,c1);
half_adder h2(s1,cin,S,c2);
assign cout=c1|c2;
endmodule

module rca_Nbit #(parameter N=16) (a, b, cin, S, cout);
input [N-1:0] a,b;
input cin;
output [N-1:0] S;
wire c[N-1:0];
output cout;

generate
    full_adder fa0(a[0],b[0],cin,S[0],c[0]);
    genvar i;
    for(i=1;i<N;i=i+1) begin
        full_adder fa1(.a(a[i]),.b(b[i]),.cin(c[i-1]),.S(S[i]),.cout(c[i]));
    end
    assign cout=c[N-1];
endgenerate
endmodule

//SUBTRACTOR
module sub #(parameter N=16) (a, b, cin, S, cout);
input [N-1:0] a,b;
wire [N-1:0] b_new;
input cin;
output [N-1:0] S;
wire c[N-1:0];
output cout;
genvar j;
for(j=0;j<N;j=j+1) begin
    assign b_new[j]=b[j]^1'b1;
end
generate
    full_adder fa0(a[0],b_new[0],cin,S[0],c[0]);
    genvar i;
    for(i=1;i<N;i=i+1) begin
        full_adder fa1(.a(a[i]),.b(b_new[i]),.cin(c[i-1]),.S(S[i]),.cout(c[i]));
    end
    assign cout=c[N-1];
endgenerate
endmodule
//FOR 2 BIT
module carry_mul(c,A,B);
input [1:0] A;
output [1:0] B;
input c;
assign B=c? A:2'b00;
endmodule
//FOR 4 BIT
module carry_mu(c,A,B);
input [3:0] A;
input c;
output [3:0] B;
assign B=c? A: 4'b0000;
endmodule
//FOR 8 BIT
module carry_m(c,A,B);
input [7:0] A;
input c;
output [7:0] B;
assign B=c? A: 8'b00000000;
endmodule

module karatsuba_2 (X, Y, Z);
input [1:0] X,Y;
output [3:0] Z;
assign Z[0]=X[0]&Y[0];
assign Z[1]=(X[1]&Y[0])^(X[0]&Y[1]);
assign Z[2]=((X[1]&Y[0])&(X[0]&Y[1]))^(X[1]&Y[1]);
assign Z[3]=((X[1]&Y[0])&(X[0]&Y[1]))&(X[1]&Y[1]);
endmodule

module karatsuba_4 (X, Y, Z);
input [3:0] X,Y;
wire [1:0] x0,x1,y0,y1,x2,y2;
assign x0=X[1:0];
assign x1=X[3:2];
assign y0=Y[1:0];
assign y1=Y[3:2];
output [7:0] Z,Z_half,Z1_half;
wire [3:0] z2,z0,z4,z3_half;
wire c1,c2;
wire c3,c4,c5,c6,c7;
wire [1:0]b1,b2;
wire [7:0] b1_ext,b2_ext,b3_ext,z3_mhalf,z3,z1,z3_fhalf,z3_half_new;
parameter N1=2;
parameter N2=4;
//z0=x0y0
karatsuba_2 k1(.X(x0),.Y(y0),.Z(z0));
//z2=x1y1
karatsuba_2 k2(.X(x1),.Y(y1),.Z(z2));
//z3=(x1+x0)(y1+y0)
rca_Nbit #(.N(2)) a2(.a(x0),.b(x1),.cin(1'b0),.S(x2),.cout(c3));
rca_Nbit #(.N(2)) a3(.a(y0),.b(y1),.cin(1'b0),.S(y2),.cout(c4));
karatsuba_2 k3(.X(x2),.Y(y2),.Z(z3_half));
carry_mul m1(.c(c3),.A(y2),.B(b1));
carry_mul m2(.c(c4),.A(x2),.B(b2));
assign b1_ext={4'b0000,b1[1:0],2'b00};
assign b2_ext={4'b0000,b2[1:0],2'b00};
assign b3_ext={3'b000,c3&c4,4'b0000};

//ADD Z3_HALF
assign z3_half_new={4'b0000,z3_half};
rca_Nbit #(.N(8)) a5(.a(z3_half_new),.b(b2_ext),.cin(1'b0),.S(z3_fhalf),.cout(c7));
rca_Nbit #(.N(8)) a4(.a(b1_ext),.b(z3_fhalf),.cin(1'b0),.S(z3_mhalf),.cout(c5));
rca_Nbit #(.N(8)) a6(.a(z3_mhalf),.b(b3_ext),.cin(1'b0),.S(z3),.cout(c6));
//z1=z3-z0-z2
wire [7:0]z0_new,z2_new;
assign z0_new={4'b0000,z0};
assign z2_new={4'b0000,z2};
wire g1,g2;
sub #(.N(8)) s1(.a(z3),.b(z0_new),.cin(1'b1),.S(Z1_half),.cout(g1));
sub #(.N(8)) s2(.a(Z1_half),.b(z2_new),.cin(1'b1),.S(z1),.cout(g2));

// end
wire [7:0] z0_e,z1_e,z2_e;
assign z0_e={4'b0,z0};
assign z1_e={z1[5:0],2'b0};
// assign z1_e=z1<<2;
assign z2_e={z2,4'b0};

rca_Nbit #(.N(8)) a0(.a(z0_e),.b(z1_e),.cin(1'b0),.S(Z_half),.cout(c1));
rca_Nbit #(.N(8)) a1(.a(Z_half),.b(z2_e),.cin(c1),.S(Z),.cout(c2));

endmodule

module karatsuba_8 (X, Y, Z);
input [7:0] X,Y;
wire [3:0] x0,x1,y0,y1,x2,y2;
assign x0=X[3:0];
assign x1=X[7:4];
assign y0=Y[3:0];
assign y1=Y[7:4];
output [15:0] Z,Z_half,Z1_half;
wire [7:0] z2,z0,z4,z3_half;
wire c1,c2;
wire c3,c4,c5,c6,c7;
wire [3:0]b1,b2;
wire [15:0] b1_ext,b2_ext,b3_ext,z3_mhalf,z3,z1,z3_fhalf,z3_half_new;
parameter N1=4;
parameter N2=8;
//z0=x0y0
karatsuba_4 k1(.X(x0),.Y(y0),.Z(z0));
//z2=x1y1
karatsuba_4 k2(.X(x1),.Y(y1),.Z(z2));
//z3=(x1+x0)(y1+y0)
rca_Nbit #(.N(4)) a2(.a(x0),.b(x1),.cin(1'b0),.S(x2),.cout(c3));
rca_Nbit #(.N(4)) a3(.a(y0),.b(y1),.cin(1'b0),.S(y2),.cout(c4));
karatsuba_4 k3(.X(x2),.Y(y2),.Z(z3_half));
carry_mu m1(.c(c3),.A(y2),.B(b1));
carry_mu m2(.c(c4),.A(x2),.B(b2));
assign b1_ext={8'b0000,b1,4'b00};
assign b2_ext={8'b0000,b2,4'b00};
assign b3_ext={7'b000,c3&c4,8'b0000};

//ADD Z3_HALF
assign z3_half_new={8'b0000,z3_half};
rca_Nbit #(.N(16)) a5(.a(z3_half_new),.b(b2_ext),.cin(1'b0),.S(z3_fhalf),.cout(c7));
rca_Nbit #(.N(16)) a4(.a(b1_ext),.b(z3_fhalf),.cin(1'b0),.S(z3_mhalf),.cout(c5));
rca_Nbit #(.N(16)) a6(.a(z3_mhalf),.b(b3_ext),.cin(1'b0),.S(z3),.cout(c6));
//z1=z3-z0-z2
wire [15:0]z0_new,z2_new;
assign z0_new={8'b0000,z0};
assign z2_new={8'b0000,z2};
wire g1,g2;
sub #(.N(16)) s1(.a(z3),.b(z0_new),.cin(1'b1),.S(Z1_half),.cout(g1));
sub #(.N(16)) s2(.a(Z1_half),.b(z2_new),.cin(1'b1),.S(z1),.cout(g2));

// end
wire [15:0] z0_e,z1_e,z2_e;
assign z0_e={8'b0,z0};
assign z1_e={z1[11:0],4'b0};
assign z2_e={z2,8'b0};

rca_Nbit #(.N(16)) a0(.a(z0_e),.b(z1_e),.cin(1'b0),.S(Z_half),.cout(c1));
rca_Nbit #(.N(16)) a1(.a(Z_half),.b(z2_e),.cin(c1),.S(Z),.cout(c2));


endmodule

module karatsuba_16 (X, Y, Z);
input [15:0] X,Y;
wire [7:0] x0,x1,y0,y1,x2,y2;
assign x0=X[7:0];
assign x1=X[15:8];
assign y0=Y[7:0];
assign y1=Y[15:8];
output [31:0] Z,Z_half,Z1_half;
wire [15:0] z2,z0,z4,z3_half;
wire c1,c2;
wire c3,c4,c5,c6,c7;
wire [7:0]b1,b2;
wire [31:0] b1_ext,b2_ext,b3_ext,z3_mhalf,z3,z1,z3_fhalf,z3_half_new;
parameter N1=8;
parameter N2=16;
//z0=x0y0
karatsuba_8 k1(.X(x0),.Y(y0),.Z(z0));  //karatsuba 1
//z2=x1y1
karatsuba_8 k2(.X(x1),.Y(y1),.Z(z2));   //karatsuba 2
//z3=(x1+x0)(y1+y0)
rca_Nbit #(.N(8)) a2(.a(x0),.b(x1),.cin(1'b0),.S(x2),.cout(c3));
rca_Nbit #(.N(8)) a3(.a(y0),.b(y1),.cin(1'b0),.S(y2),.cout(c4));

karatsuba_8 k3(.X(x2),.Y(y2),.Z(z3_half)); //karatsuba 3
carry_m m1(.c(c3),.A(y2),.B(b1));
carry_m m2(.c(c4),.A(x2),.B(b2));
assign b1_ext={16'b0000000000000000,b1,8'b00000000};
assign b2_ext={16'b0000000000000000,b2,8'b00000000};
assign b3_ext={15'b000000000000000,c3&c4,16'b0000000000000000};

//ADD Z3_HALF
assign z3_half_new={16'b0000000000000000,z3_half};
rca_Nbit #(.N(32)) a5(.a(z3_half_new),.b(b2_ext),.cin(1'b0),.S(z3_fhalf),.cout(c7));
rca_Nbit #(.N(32)) a4(.a(b1_ext),.b(z3_fhalf),.cin(1'b0),.S(z3_mhalf),.cout(c5));
rca_Nbit #(.N(32)) a6(.a(z3_mhalf),.b(b3_ext),.cin(1'b0),.S(z3),.cout(c6));
//z1=z3-z0-z2
wire [31:0]z0_new,z2_new;
assign z0_new={16'b0000000000000000,z0};
assign z2_new={16'b0000000000000000,z2};
wire g1,g2;
sub #(.N(32)) s1(.a(z3),.b(z0_new),.cin(1'b1),.S(Z1_half),.cout(g1));
sub #(.N(32)) s2(.a(Z1_half),.b(z2_new),.cin(1'b1),.S(z1),.cout(g2));

// end
wire [31:0] z0_e,z1_e,z2_e;
assign z0_e={16'b0,z0};
assign z1_e={z1[23:0],8'b0};
assign z2_e={z2,16'b0};

rca_Nbit #(.N(32)) a0(.a(z0_e),.b(z1_e),.cin(1'b0),.S(Z_half),.cout(c1));
rca_Nbit #(.N(32)) a1(.a(Z_half),.b(z2_e),.cin(c1),.S(Z),.cout(c2));

endmodule