import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { GUI } from 'dat.gui';

// Define constants
const DEFAULT_BG_COLOR = '#888888';
const DEFAULT_FONT_COLOR = '#ffffff';
const DEFAULT_FONT_SIZE = 14;
const MIN_FONT_SIZE = 8;
const MAX_FONT_SIZE = 30;
const DEFAULT_RUB_RADIUS = 8;
const MIN_RUB_RADIUS = 1;
const MAX_RUB_RADIUS = 45;
const DEFAULT_OPACITY_LEVEL = 0.3;

// Initialize state variables
let scene, camera, renderer, controls, gui;
let rubs, rubShapeData = [], geometry, materials = [];
let params;


function updateTagPosition(rub, tagElement) {
  if (tagElement === null)
    return;
  
  const vector = new THREE.Vector3();
  rub.getWorldPosition(vector);
  vector.project(camera);
  const x = (vector.x * 0.5 + 0.5) * window.innerWidth;
  const y = (-vector.y * 0.5 + 0.5) * window.innerHeight;
  tagElement.style.left = `${x}px`;
  tagElement.style.top = `${y}px`;
  tagElement.style.display = (rub.visible && params.showTags) ? 'block' : 'none';
}

// Update function to apply font color to all tags
function updateFontColor() {
  rubShapeData.forEach(({ tagElement }) => {
    if (tagElement) {
      tagElement.style.color = params.fontColor;
    }
  });
}

// Update function to apply font size to all tags
function updateFontSize() {
  rubShapeData.forEach(({ tagElement }) => {
    if (tagElement) {
      tagElement.style.fontSize = `${params.fontSize}px`;
    }
  });
}

function createTagElement(tag) {
  const tagElement = document.createElement('div');
  tagElement.style.position = 'absolute';
  tagElement.style.color = params.fontColor;
  tagElement.textContent = tag;
  document.body.appendChild(tagElement);
  return tagElement;
}

function createRubisco(r, material, centerOfMass) {
  const rub = new THREE.Mesh(geometry, material);

  const tdrot = THREE.MathUtils.degToRad(r.tdrot);
  const tilt = THREE.MathUtils.degToRad(r.tilt);
  const narot = THREE.MathUtils.degToRad(r.narot);

  rub.rotateZ(-tdrot).rotateX(-tilt).rotateZ(-narot);
  rub.rotateX(THREE.MathUtils.degToRad(90));

  rub.position.set(r.x - centerOfMass[0], r.y - centerOfMass[1], r.z - centerOfMass[2]);
  return rub;
}

function generateRandomColor() {
  return Math.floor(Math.random() * 16777215);
}

function calculateCenterOfMass(rubs) {
  const centerOfMass = [0, 0, 0];
  rubs.forEach(r => {
    centerOfMass[0] += r.x;
    centerOfMass[1] += r.y;
    centerOfMass[2] += r.z;
  });
  return centerOfMass.map(c => c / rubs.length);
}

function drawCarb(rubs, tagsOfInterest) {
  // Clear existing objects from scene
  rubShapeData.forEach(({ rub, tagElement }) => {
    scene.remove(rub);
    tagElement && document.body.removeChild(tagElement);
  });
  rubShapeData = [];
  
  // Boolean indicator of whether user specified multiple chains or individual Rubisco to highlight
  const multiple = tagsOfInterest && tagsOfInterest.every(Array.isArray);
  if (multiple && params.colorMode === 'Random') {
    params.randomChainColors = tagsOfInterest.map(() => generateRandomColor());
  }

  geometry = new THREE.CylinderGeometry(params.radius, params.radius, 45, 4, 1);

  rubs.forEach(r => {
    let color;
    let displayTag = false;
    if (tagsOfInterest === null) {
      color = params.colorMode === 'Random' ? generateRandomColor() : params.uniformColor;
      displayTag = true;
    } else if (multiple) {
      const chainIndex = tagsOfInterest.findIndex(chain => chain.includes(r.tag));
      color = chainIndex !== -1 ? (params.colorMode === 'Random' ? params.randomChainColors[chainIndex] : params.uniformColor) : '#000000';
      displayTag = chainIndex !== -1;
    } else {
      color = tagsOfInterest.includes(r.tag) ? (params.colorMode == 'Random' ? generateRandomColor() : params.uniformColor) : '#000000';
      displayTag = tagsOfInterest.includes(r.tag);
    }

    const material = new THREE.MeshBasicMaterial({
      color,
      transparent: true,
      opacity: displayTag ? 1 : params.opacity,
    });

    const centerOfMass = calculateCenterOfMass(rubs);
    const rub = createRubisco(r, material, centerOfMass);
    rub.userData.tag = r.tag;
    rub.visible = !(params.hiddenTags.includes(r.tag));
    const tagElement = displayTag ? createTagElement(r.tag) : null;
    rubShapeData.push({ rub, tagElement, updateTagPosition: () => {
      if (tagElement) {
        updateTagPosition(rub, tagElement);
      }
    }});
    scene.add(rub);
  });
}

function updateRadius() {
  rubShapeData.forEach(({ rub }) => {
    rub.geometry.dispose(); // Dispose of old geometry to prevent memory leaks
    rub.geometry = new THREE.CylinderGeometry(params.radius, params.radius, 45, 4, 1);
  });
}

function updateOpacity() {
  rubShapeData.forEach(({ rub }) => {
    if (rub.material.opacity !== 1) {
      rub.material.opacity = params.opacity;
    }
  });
}

function createOverlayText() {
  const overlayText = document.createElement('div');
  overlayText.style.position = 'absolute';
  overlayText.style.top = '16px'; // Adjust the position
  overlayText.style.left = '15px'; // Adjust the position
  overlayText.style.backgroundColor = 'rgba(0, 0, 0, 0.7)'; // Semi-transparent background
  overlayText.style.color = '#ffffff'; // Text color
  overlayText.style.padding = '10px';
  overlayText.style.fontSize = '14px'; // Adjust font size
  overlayText.style.zIndex = '100'; // Ensure it appears above other elements
  overlayText.innerHTML = `
    <strong>Keyboard Controls:</strong><br>
    <span style="font-size: 12px;">
      <strong>'R':</strong> Re-generate random colors<br>
      <strong>'T':</strong> Toggle tag visibility
    </span>
  `;
  document.body.appendChild(overlayText);
}

// Function to create overlay control for showing/hiding Rubisco objects
function createVisibilityControl() {
  const controlBox = document.createElement('div');
  controlBox.style.position = 'absolute';
  controlBox.style.bottom = '16px'; // Position at the bottom right
  controlBox.style.right = '15px';
  controlBox.style.backgroundColor = 'rgba(0, 0, 0, 0.7)';
  controlBox.style.color = '#ffffff';
  controlBox.style.padding = '10px';
  controlBox.style.fontSize = '14px';
  controlBox.style.maxWidth = '350px';
  controlBox.style.overflowWrap = 'break-word';

  // Create input for tag entry
  const tagInput = document.createElement('input');
  tagInput.type = 'number';
  tagInput.placeholder = 'Enter a tag';

  // Create button for toggling visibility
  const toggleButton = document.createElement('button');
  toggleButton.textContent = 'Hide/Show';
  toggleButton.style.marginLeft = '4px';

  // Create button to reset to all rubisco visible
  const resetButton = document.createElement('button');
  resetButton.textContent = 'Reset';
  resetButton.style.marginLeft = '4px';

  // Create list to display hidden tags
  const hiddenTagsDisplay = document.createElement('div');
  hiddenTagsDisplay.textContent = "Hidden Tags: " + params.hiddenTags;
  hiddenTagsDisplay.style.marginTop = '10px';

  // Append elements to the control box
  controlBox.appendChild(tagInput);
  controlBox.appendChild(toggleButton);
  controlBox.appendChild(resetButton);
  controlBox.appendChild(hiddenTagsDisplay);
  document.body.appendChild(controlBox);

  // Toggle visibility function
  function toggleTagVisibility() {
    const tag = parseInt(tagInput.value);

    const itemIndex = rubShapeData.findIndex(item => item.rub.userData.tag === tag);

    if (itemIndex !== -1) {
      const { rub } = rubShapeData[itemIndex];

      if (params.hiddenTags.includes(tag.toString())) {
        // If the tag is hidden, show it
        rub.visible = true;
        params.hiddenTags = params.hiddenTags.split(', ').filter(t => t !== tag.toString()).join(', ');
      } else {
        // If the tag is visible, hide it
        rub.visible = false;
	params.hiddenTags += params.hiddenTags.length > 0 ? `, ${tag}` : `${tag}`;
      }

      // Update the hidden tags display string
      hiddenTagsDisplay.textContent = "Hidden Tags: " + params.hiddenTags;
    } else {
      alert("Tag not found or invalid.");
    }

    // Clear the input field
    tagInput.value = '';
  }

  // Add click event for toggle button
  toggleButton.addEventListener('click', toggleTagVisibility);

  // Add keydown event for Enter key on the input field
  tagInput.addEventListener('keydown', (event) => {
    if (event.key === 'Enter') {
      toggleTagVisibility();
    }
  });

  // Add click event for reset button
  resetButton.addEventListener('click', () => {
    rubShapeData.forEach(({ rub }) => {
      rub.visible = true;
    });
    params.hiddenTags = "";
    hiddenTagsDisplay.textContent = "Hidden Tags: ";
  });
}

function animate() {
  requestAnimationFrame(animate);
  rubShapeData.forEach(({ updateTagPosition }) => updateTagPosition());
  renderer.render(scene, camera);
}

async function main() {
	const filename = prompt("Please enter the name of the .tbl file to read: ");

	const header = ["tag", "aligned", "averaged", "dx", "dy", "dz", "tdrot", "tilt", "narot", "cc", "cc2", "cpu", "ftype", "ymintilt", "ymaxtilt", "xmintilt", "xmaxtilt", "fs1", "fs2", "tomo", "reg", "class", "annotation", "x", "y", "z", "dshift", "daxis", "dnarot", "dcc", "otag", "npar", "ref", "sref", "apix", "def", "eig1", "eig2"];

	// Read input data from .tbl file
	const response = await fetch(filename);
	const data = await response.text();
	const lines = data.trim().split('\n');
	const allData = lines.map(line => {
	  const fields = line.split(' ');
	  const obj = {};
	  header.forEach((h, i) => obj[h] = parseFloat(fields[i]));
	  return obj;
	}); 

	// Ask for Carboxysome Number to model
	const input = prompt("Enter the Carboxysome no. that you wish to model: ");
	const carbIndex = parseInt(input);

	// Extract Rubisco in the specified Carboxysome
	rubs = allData.filter(rub => rub.reg === carbIndex);
	console.log(rubs);

	// Ask for tags of Rubisco to highlight
	const tags = prompt("Enter the tags of the Rubisco that you wish to highlight. You may enter a list of tags (e.g., [1,2,3]) or a list of lists of tags to view them as chains (e.g., [[1,2],[2,3]]). Otherwise, click cancel to color all rubisco with random colors.");
	let tagsOfInterest = null;
	if (tags !== "") {
	  tagsOfInterest = JSON.parse(tags);
	}

	// Set up 3D Scene
	scene = new THREE.Scene();
	scene.background = new THREE.Color(DEFAULT_BG_COLOR);

	camera = new THREE.OrthographicCamera(window.innerWidth / -2, window.innerWidth / 2, window.innerHeight / 2, window.innerHeight / - 2, 0.1, 1000);
	camera.position.z = 400;
	

	// Initialize modifiable parameters
	params = {
	  backgroundColor: DEFAULT_BG_COLOR,
	  fontColor: DEFAULT_FONT_COLOR,
	  fontSize: DEFAULT_FONT_SIZE,
	  radius: DEFAULT_RUB_RADIUS,
	  opacity: DEFAULT_OPACITY_LEVEL,
	  colorMode: 'Random',
	  randomChainColors: [],
	  uniformColor: '#ff0000',
	  hiddenTags: "",
	};


	// Initialize User Controls GUI
	gui = new GUI({width: 500});
	// Add background color picker
	gui.addColor(params, 'backgroundColor').onChange(value => scene.background.set(value));
	gui.addColor(params, 'fontColor').name('Tag Font Color').onChange(updateFontColor);
	gui.add(params, 'fontSize', MIN_FONT_SIZE, MAX_FONT_SIZE).name('Tag Font Size').onChange(updateFontSize);
	// Add slider for radius of Rubisco objects
	gui.add(params, 'radius', MIN_RUB_RADIUS, MAX_RUB_RADIUS).onChange(updateRadius);

	// If not all rubisco are highlighted, add a slider for the opacity of non-marked rubisco
	if (tagsOfInterest !== null) {
	  gui.add(params, 'opacity', 0, 0.99).onChange(updateOpacity);
	}

	const colorFolder = gui.addFolder('Rubisco Color Options');
	colorFolder.open();
	colorFolder.add(params, 'colorMode', ['Random', 'Choose Color']).name('Color Mode').onChange((value) => {
    	  if (value === 'Choose Color') {
      	    // Show uniform color picker if "Choose Color" mode is selected
      	    uniformColorController = colorFolder.addColor(params, 'uniformColor').name('Uniform Color').onChange(() => {drawCarb(rubs, tagsOfInterest);});
    	  } else if (uniformColorController) {
      	    // Hide the color picker if switching back to "Random"
      	    colorFolder.remove(uniformColorController);
      	    uniformColorController = null;
    	  }
    	  drawCarb(rubs, tagsOfInterest);  // Update colors when mode changes
  	});
	
	let uniformColorController = null;
	
	// Set up renderer
	renderer = new THREE.WebGLRenderer();
	renderer.setSize(window.innerWidth, window.innerHeight);
	document.body.appendChild(renderer.domElement);

	// Set up Orbit Controls (allows the user to drag, move, zoom in the model)
	controls = new OrbitControls(camera, renderer.domElement);
	
	
	// Draw Carboxysome
	drawCarb(rubs, tagsOfInterest);
	
	// Start animation
	animate();

	// Create Visibility controls
	createVisibilityControl();

	// Set up keyboard controls
	createOverlayText();
	document.addEventListener('keydown', event => {
  	  const key = event.key.toLowerCase();
  	  if (key === 'r' && params.colorMode !== 'Choose Color') {
	    drawCarb(rubs, tagsOfInterest);
  	  }
  	  if (key === 't') {
	    params.showTags = !params.showTags;
  	  }
	});
}

// Run the main function to initiate the program
main();
