import processing.pdf.*;

void setup() {
  size(500, 800, PDF, "flowchart-orfansingmult.pdf");
  background(255);
}

void polygon(int n, float cx, float cy, float r) {
  float angle = 360.0 / n;

  beginShape();
  for (int i = 0; i < n; i++) {
    vertex(cx + r * cos(radians(angle * i)),
      cy + r * sin(radians(angle * i)));
  }
  endShape(CLOSE);
}

void arrow(int x1, int y1, int x2, int y2) {
  line(x1, y1, x2, y2);
  pushMatrix();
  translate(x2, y2);
  float a = atan2(x1-x2, y2-y1);
  rotate(a);
  line(0, 0, -10, -10);
  line(0, 0, 10, -10);
  popMatrix();
} 


void draw() {
  fill(255);
  polygon(4, 380, 120, 120);
  fill(0);
  textSize(21);
  text("Is the number",310,90);
  text("of species in the",310,110);
  text("cluster less than",310,130);
  text("or equal to 3?",310,150);
  
  arrow(380, 240, 380, 350);
  pushMatrix();
  translate(380,280);
  rotate(HALF_PI);
  text("No", 0, 0);
  popMatrix();
  
  arrow(260, 120, 150, 120);
  text("Yes",200 , 120);
  
  fill(255);
  polygon(4, 380, 470, 120);
  fill(0);
  textSize(21);
  text("Are less than",310,450);
  text("30% of the genes",310,470);
  text("in the cluster",310,490);
  text("duplicated?",310,510);
  
  arrow(380, 590, 380, 700);
  pushMatrix();
  translate(380,630);
  rotate(HALF_PI);
  text("No", 0, 0);
  popMatrix();
  
  arrow(260, 470, 150, 470);
  text("Yes",200 , 470);
  
  fill(#3FE31B);
  rect(0,80,150,80);
  fill(0);
  text("ORFan", 45, 125);
  
  fill(#DEA324);
  rect(0,430,150,80);
  fill(0);
  text("Single copy", 15, 475);
  
  fill(#9B9B9B);
  rect(305,700,150,80);
  fill(0);
  text("Multiple copy", 310, 745);
  
  exit();
}