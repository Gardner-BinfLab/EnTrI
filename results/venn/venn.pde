import processing.pdf.*;

void setup() {
  size(560, 400, PDF, "venn.pdf");
}

void draw() {
  fill(200, 50, 100, 80);
  ellipse(200, 200, 400, 400);
  
  fill(100, 50, 200, 80);
  ellipse(360, 200, 400, 400);
  
  fill(200, 100, 50, 80);
  ellipse(200, 200, 240, 240);
  
  fill(50, 100, 200, 80);
  ellipse(360, 200, 240, 240);
  
  fill(0);
  textSize(40);

  pushMatrix();
  translate(50, 120);
  rotate(HALF_PI);
  text("Pres-Abs", 0, 0);
  popMatrix();

  pushMatrix();
  translate(500, 120);
  rotate(HALF_PI);
  text("Abs-Pres", 0, 0);
  popMatrix();

  textSize(30);
  text("Pres-", 250, 60);
  text("Pres", 250, 90);
  textSize(40);

  pushMatrix();
  translate(120, 120);
  rotate(HALF_PI);
  text("Ess-Abs", 0, 0);
  popMatrix();

  pushMatrix();
  translate(410, 120);
  rotate(HALF_PI);
  text("Abs-Ess", 0, 0);
  popMatrix();
  
  pushMatrix();
  translate(270, 130);
  rotate(HALF_PI);
  text("Ess-Ess", 0, 0);
  popMatrix();
  
  pushMatrix();
  translate(190, 120);
  rotate(HALF_PI);
  text("Ess-NoEs", 0, 0);
  popMatrix();
  
  pushMatrix();
  translate(330, 120);
  rotate(HALF_PI);
  text("NoEs-Ess", 0, 0);
  popMatrix();
  
  exit();
}