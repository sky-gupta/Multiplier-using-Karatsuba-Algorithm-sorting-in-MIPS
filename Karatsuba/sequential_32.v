/* 32-bit simple karatsuba multiplier */

/*32-bit Karatsuba multipliction using a single 16-bit module*/

module iterative_karatsuba_32_16(clk, rst, enable, A, B, C);
    input clk;
    input rst;
    input [31:0] A;
    input [31:0] B;
    output [63:0] C;
    
    input enable;
    
    
    wire [1:0] sel_x;
    wire [1:0] sel_y;
    
    wire [1:0] sel_z;
    wire [1:0] sel_T;
    
    
    wire done;
    wire en_z;
    wire en_T;
    
    
    wire [32:0] h1;
    wire [32:0] h2;
    wire [63:0] g1;
    wire [63:0] g2;
    
    assign C = g2;
    reg_with_enable #(.N(63)) Z(.clk(clk), .rst(rst), .en(en_z), .X(g1), .O(g2) );  // Fill in the proper size of the register
    reg_with_enable #(.N(32)) T(.clk(clk), .rst(rst), .en(en_T), .X(h1), .O(h2) );  // Fill in the proper size of the register
    
    iterative_karatsuba_datapath dp(.clk(clk), .rst(rst), .X(A), .Y(B), .Z(g2), .T(h2), .sel_x(sel_x), .sel_y(sel_y), .sel_z(sel_z), .sel_T(sel_T), .en_z(en_z), .en_T(en_T), .done(done), .W1(g1), .W2(h1));
    iterative_karatsuba_control control(.clk(clk),.rst(rst), .enable(enable), .sel_x(sel_x), .sel_y(sel_y), .sel_z(sel_z), .sel_T(sel_T), .en_z(en_z), .en_T(en_T), .done(done));
    
endmodule

module iterative_karatsuba_datapath(clk, rst, X, Y, T, Z, sel_x, sel_y, en_z, sel_z, en_T, sel_T, done, W1, W2);
    input clk;
    input rst;
    input [31:0] X;    // input X
    input [31:0] Y;    // Input Y
    input [32:0] T;    // input which sums X_h*Y_h and X_l*Y_l (its also a feedback through the register)
    input [63:0] Z;    // input which calculates the final outcome (its also a feedback through the register)
    output [63:0] W1;  // Signals going to the registers as input
    output [32:0] W2;  // signals hoing to the registers as input
    

    input [1:0] sel_x;  // control signal 
    input [1:0] sel_y;  // control signal 
    
    input en_z;         // control signal 
    input [1:0] sel_z;  // control signal 
    input en_T;         // control signal 
    input [1:0] sel_T;  // control signal 
    
    input done;         // Final done signal
    
    
   wire [15:0] x0,x1,y0,y1,x2,y2;
    assign x0=X[15:0];
    assign x1=X[31:16];
    assign y0=Y[15:0];
    assign y1=Y[31:16];
    wire [63:0] Z_final,Z_half,Z1_half;
    wire [31:0] z2,z0,z4,z3_half;
    wire c1,c2;
    wire c3,c4,c5,c6,c7;
    wire [15:0]b1,b2;
    wire [63:0] b1_ext,b2_ext,b3_ext,z3_mhalf,z3,z1,z3_fhalf,z3_half_new;
    parameter N1=16;
    parameter N2=32;

    wire [15:0]copy_x,copy_y;
    wire [31:0]copy_z;
    mult_16 m(.X(copy_x),.Y(copy_y),.Z(copy_z));
    assign copy_x=(sel_x==2'b00)? x0:((sel_x==2'b01)? x1: x2);
    assign copy_y=(sel_x==2'b00)? y0:((sel_x==2'b01)? y1: y2);
    assign z0=(sel_x==2'b00)? copy_z:z0;
    assign z2=(sel_x==2'b01)? copy_z:z2;
    assign z3_half=(sel_x==2'b10)? copy_z:z3_half;
    //z3=(x1+x0)(y1+y0)
    rca_Nbit #(.N(16)) a2(.a(x0),.b(x1),.cin(1'b0),.S(x2),.cout(c3));
    rca_Nbit #(.N(16)) a3(.a(y0),.b(y1),.cin(1'b0),.S(y2),.cout(c4));
    carry_m m1(.c(c3),.A(y2),.B(b1));
    carry_m m2(.c(c4),.A(x2),.B(b2));
    assign b1_ext={32'b0000000000000000,b1,16'b00000000};
    assign b2_ext={32'b0000000000000000,b2,16'b00000000};
    assign b3_ext={31'b000000000000000,c3&c4,32'b0000000000000000};

    //ADD Z3_HALF
    assign z3_half_new={16'b0000000000000000,z3_half};
    rca_Nbit #(.N(64)) a5(.a(z3_half_new),.b(b2_ext),.cin(1'b0),.S(z3_fhalf),.cout(c7));
    rca_Nbit #(.N(64)) a4(.a(b1_ext),.b(z3_fhalf),.cin(1'b0),.S(z3_mhalf),.cout(c5));
    rca_Nbit #(.N(64)) a6(.a(z3_mhalf),.b(b3_ext),.cin(1'b0),.S(z3),.cout(c6));
    //z1=z3-z0-z2
    wire [63:0]z0_new,z2_new;
    assign z0_new={32'b0000000000000000,z0};
    assign z2_new={32'b0000000000000000,z2};
    wire g1,g2;
    sub #(.N(64)) s1(.a(z3),.b(z0_new),.cin(1'b1),.S(Z1_half),.cout(g1));
    sub #(.N(64)) s2(.a(Z1_half),.b(z2_new),.cin(1'b1),.S(z1),.cout(g2));

    // end
    wire [63:0] z0_e,z1_e,z2_e;
    assign z0_e={32'b0,z0};
    assign z1_e={z1[47:0],16'b0};
    assign z2_e={z2,32'b0};

    rca_Nbit #(.N(64)) a0(.a(z0_e),.b(z1_e),.cin(1'b0),.S(Z_half),.cout(c1));
    rca_Nbit #(.N(64)) a1(.a(Z_half),.b(z2_e),.cin(c1),.S(Z_final),.cout(c2));
    assign W1=Z_final;
    //-------------------------------------------------------------------------------------------------
    //--------------------------------------------------------

endmodule


module iterative_karatsuba_control(clk,rst, enable, sel_x, sel_y, sel_z, sel_T, en_z, en_T, done);
    input clk;
    input rst;
    input enable;
    
    output reg [1:0] sel_x;
    output reg [1:0] sel_y;
    
    output reg [1:0] sel_z;
    output reg [1:0] sel_T;    
    
    output reg en_z;
    output reg en_T;
    
    
    output reg done;
    
    reg [5:0] state, nxt_state;
    parameter S0 = 6'b000001;   // initial state
   // <define the rest of the states here>
    parameter S1=6'b000010;
    parameter S2=6'b000100;
    parameter S3=6'b001000;

    always @(posedge clk) begin
        if (rst) begin
            state <= S0;
        end
        else if (enable) begin
            state <= nxt_state;
        end
    end
    

    always@(*) begin
        case(state) 
            S0:
                begin
                    sel_x=2'b00;
                    sel_y=2'b00;
                    sel_T=2'b00;
                    nxt_state=S1;
                    en_z=1;
                end
            S1: 
                begin
					// Write your output and next state equations here
                    sel_x=2'b01;
                    sel_y=2'b01;
                    sel_T=2'b01;
                    nxt_state=S2;
                    en_z=1;
                end
            S2:
                begin
                    sel_x=2'b10;
                    sel_y=2'b10;
                    sel_T=2'b10;
                    en_z=1;
                    nxt_state=S3;
                end
            S3:
                begin
                    done=1;
                    en_z=1;
                    nxt_state=S3;
                end
			// Define the rest of the states
            default: 
				// Don't forget the default   
                begin   
                sel_x=2'b00;
                en_z=1;
                nxt_state=S0;  
                end
        endcase
        
    end

endmodule


module reg_with_enable #(parameter N = 32) (clk, rst, en, X, O );
    input [N:0] X;
    input clk;
    input rst;
    input en;
    output [N:0] O;
    
    reg [N:0] R;
    
    always@(posedge clk) begin
        if (rst) begin
            R <= {N{1'b0}};
        end
        if (en) begin
            R <= X;
        end
    end
    assign O = R;
endmodule


//My helpers
module carry_m(c,A,B);
input [15:0] A;
input c;
output [15:0] B;
assign B=c? A: 16'b00000000;
endmodule
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


/*-------------------Supporting Modules--------------------*/
/*------------- Iterative Karatsuba: 32-bit Karatsuba using a single 16-bit Module*/

module mult_16(X, Y, Z);
input [15:0] X;
input [15:0] Y;
output [31:0] Z;

assign Z = X*Y;

endmodule


module check_subtract (A, B, C);
 input [7:0] A;
 input [7:0] B;
 output [8:0] C;
 
 assign C = A - B; 
endmodule



/* N-bit RCA adder (Unsigned) */
module adder_Nbit #(parameter N = 32) (a, b, cin, S, cout);
input [N-1:0] a;
input [N-1:0] b;
input cin;
output [N-1:0] S;
output cout;

wire [N:0] cr;  

assign cr[0] = cin;


generate
    genvar i;
    for (i = 0; i < N; i = i + 1) begin
        full_adder addi (.a(a[i]), .b(b[i]), .cin(cr[i]), .S(S[i]), .cout(cr[i+1]));
    end
endgenerate    


assign cout = cr[N];

endmodule


module Not_Nbit #(parameter N = 32) (a,c);
input [N-1:0] a;
output [N-1:0] c;

generate
genvar i;
for (i = 0; i < N; i = i+1) begin
    assign c[i] = ~a[i];
end
endgenerate 

endmodule


/* 2's Complement (N-bit) */
module Complement2_Nbit #(parameter N = 32) (a, c, cout_comp);

input [N-1:0] a;
output [N-1:0] c;
output cout_comp;

wire [N-1:0] b;
wire ccomp;

Not_Nbit #(.N(N)) compl(.a(a),.c(b));
adder_Nbit #(.N(N)) addc(.a(b), .b({ {N-1{1'b0}} ,1'b1 }), .cin(1'b0), .S(c), .cout(ccomp));

assign cout_comp = ccomp;

endmodule


/* N-bit Subtract (Unsigned) */
module subtract_Nbit #(parameter N = 32) (a, b, cin, S, ov, cout_sub);

input [N-1:0] a;
input [N-1:0] b;
input cin;
output [N-1:0] S;
output ov;
output cout_sub;

wire [N-1:0] minusb;
wire cout;
wire ccomp;

Complement2_Nbit #(.N(N)) compl(.a(b),.c(minusb), .cout_comp(ccomp));
adder_Nbit #(.N(N)) addc(.a(a), .b(minusb), .cin(1'b0), .S(S), .cout(cout));

assign ov = (~(a[N-1] ^ minusb[N-1])) & (a[N-1] ^ S[N-1]);
assign cout_sub = cout | ccomp;

endmodule



/* n-bit Left-shift */

module Left_barrel_Nbit #(parameter N = 32)(a, n, c);

input [N-1:0] a;
input [$clog2(N)-1:0] n;
output [N-1:0] c;


generate
genvar i;
for (i = 0; i < $clog2(N); i = i + 1 ) begin: stage
    localparam integer t = 2**i;
    wire [N-1:0] si;
    if (i == 0) 
    begin 
        assign si = n[i]? {a[N-t:0], {t{1'b0}}} : a;
    end    
    else begin 
        assign si = n[i]? {stage[i-1].si[N-t:0], {t{1'b0}}} : stage[i-1].si;
    end
end
endgenerate

assign c = stage[$clog2(N)-1].si;

endmodule



