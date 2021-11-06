const Settings_ElementTypeBit = 0;
const Settings_ElementTypeMask = (1 << 4) - 1;
const Settings_ConstructBlockFromSettings = (1 << 4);
const Settings_ForceResetBlock = (1 << 5);
const Settings_EnergyBit = 6;
const Settings_EnergyMask = (1 << 5) - 1;
const Settings_XpbdSolveBit = 11;
const Settings_XpbdSolveMask = 1;
const Settings_ShapeBit = 12;
const Settings_ShapeMask = (1 << 4) - 1;
const Settings_ElementPatternBit = 16;
const Settings_ElementPatternMask = (1 << 2) - 1;
const Settings_ConstraintOrderBit = 18;
const Settings_ConstraintOrderMask = (1 << 2) - 1;
const Settings_RayleighTypeBit = 20;
const Settings_RayleighTypeMask = (1 << 2) - 1;
const Settings_Rotate90Degrees = (1 << 24);
const Settings_Centered = (1 << 25);
const Settings_LockLeft = (1 << 26);
const Settings_LockRight = (1 << 27);
const Settings_RotateLock = (1 << 28);

const Element_Null = 0;
const Element_T3 = 1;
const Element_T6 = 2;
const Element_Q4 = 3;
const Element_Q9 = 4;
const Element_T4 = 5;
const Element_T10 = 6;
const Element_H8 = 7;
const Element_H27 = 8;

const Energy_Pixar = 0;
const Energy_PixarReduced = 1;
const Energy_PixarSel = 2;
const Energy_Mixed = 3;
const Energy_MixedSel = 4;
const Energy_YeohSkin = 5;
const Energy_YeohSkinSel = 6;
const Energy_YeohSkinFast = 7;
const Energy_ContinuousPixar = 8;
const Energy_ContinuousMixed = 9;
const Energy_ContinuousSkin = 10;
const Energy_CubeNeo = 11;
const Energy_CubeSkin = 12;

const Shape_Single = 0;
const Shape_Line = 1;
const Shape_BeamL = 2;
const Shape_BeamM = 3;
const Shape_BeamH = 4;
const Shape_BeamL1x2 = 5;
const Shape_BeamL2x1 = 6;
const Shape_BeamL4x1 = 7;
const Shape_BeamL8x1 = 8;
const Shape_BoxL = 9;
const Shape_BoxM = 10;
const Shape_BoxH = 11;
const Shape_Armadillo = 12;

const Pattern_Uniform = 0;
const Pattern_Mirrored = 1;

const Rayleigh_Paper = 0;
const Rayleigh_Limit = 1;
const Rayleigh_Post = 2;
const Rayleigh_PostAmortized = 3;

function makeDefaultSettings(element) {
	return {
	timeScale: 1.0,
	stepsPerSecond: 3000.0,
	volumePasses: 0,
	gravity: 0.05,
	compliance: 1.0,
	damping: 0.005,
	pbdDamping: 0.03,
	drag: 0.002,
	poissonsRatio: 0.495,
	leftRightSeparation: 1.0,
	wonkiness: 0.0,
	flags: Settings_ConstructBlockFromSettings |
		(element << Settings_ElementTypeBit) |
		(Energy_MixedSel << Settings_EnergyBit) |
		(1 << Settings_XpbdSolveBit) |
		(Shape_BoxL << Settings_ShapeBit) |
		(Pattern_Uniform << Settings_ElementPatternBit) |
		(Rayleigh_PostAmortized << Settings_RayleighTypeBit) |
		Settings_LockLeft,
	};
}
var settings = makeDefaultSettings(Element_Q4);
var settings1 = makeDefaultSettings(Element_Null);
const SettingsId = 'settings_v4';
const settings0Json = localStorage.getItem(`${SettingsId}[0]`);
if (settings0Json) { settings = JSON.parse(settings0Json); }
const settings1Json = localStorage.getItem(`${SettingsId}[1]`);
if (settings1Json) { settings1 = JSON.parse(settings1Json); }
var links = [];

function updateSettings(simIdx) {
	simIdx = simIdx ? 1 : 0;
	let pushSettings = simIdx ? settings1 : settings;
	rpc('UpdateSettings',
		simIdx ? 1 : 0,
		pushSettings.timeScale,
		pushSettings.stepsPerSecond,
		pushSettings.volumePasses,
		pushSettings.gravity,
		pushSettings.compliance,
		pushSettings.damping,
		pushSettings.pbdDamping,
		pushSettings.drag,
		pushSettings.poissonsRatio,
		pushSettings.wonkiness,
		pushSettings.leftRightSeparation,
		pushSettings.flags
	);
	window.localStorage.setItem(`${SettingsId}[${simIdx}]`, JSON.stringify(pushSettings));
}

function resetBlocks() {
	settings.flags |= Settings_ForceResetBlock;
	settings1.flags |= Settings_ForceResetBlock;
	updateSettings(0);
	updateSettings(1);
	settings.flags &= ~Settings_ForceResetBlock;
	settings1.flags &= ~Settings_ForceResetBlock;
}
function resetSettings() {
	settings = makeDefaultSettings(Element_Q4);
	settings1 = makeDefaultSettings(Element_Null);
	resetBlocks();
	location.reload();
	return false;
}

function toggleAllLinks() {
	let target = !links[3].linkCheckbox.checked;
	for (let i = 3; i < links.length; i++) {
		links[i].linkCheckbox.checked = target;
		links[i].linkCheckbox.dispatchEvent(new Event('change'));
	}
}

function addElement(parentElement, type, className) {
	let element = document.createElement(type);
	if (className) {
		element.className = className;
	}
	parentElement.appendChild(element);
	return element;
}
function addLabel(parentElement, className, text) {
	let label = addElement(parentElement, 'label', className);
	label.innerHTML = text;
	return label;
}
function addInput(parentElement, className, type) {
	let input = addElement(parentElement, 'input', className);
	input.type = type;
	return input;
}

class CheckboxControl {
	constructor(parentElement, flags, labels) {
		console.assert(flags.length == labels.length);
		this.flags = flags;
		this.onchange = null;

		this.rootElement = addElement(parentElement, 'div', 'control-input-bag');
		let row = addElement(this.rootElement, 'div', 'control-input-bag-row');
		this.checkboxes = [];
		for (let i = 0; i < flags.length; i++) {
			let labelNode = addLabel(row, 'control-input-bag-contents', labels[i]);
			let checkbox = addInput(labelNode, '', 'checkbox');
			checkbox.addEventListener('change', () => {
				if (this.onchange) { this.onchange(checkbox.checked); }
			});
			this.checkboxes[i] = checkbox;
		}
	}
	storeState() {
		let state = [];
		for (let i = 0; i < this.checkboxes.length; i++) {
			state[i] = this.checkboxes[i].checked;
		}
		return state;
	}
	loadState(state) {
		for (let i = 0; i < this.checkboxes.length; i++) {
			this.checkboxes[i].checked = state[i];
		}
	}
	storeToSettings(settings) {
		let changed = false;
		for (let i = 0; i < this.checkboxes.length; i++) {
			if (this.checkboxes[i].checked != !!(settings.flags & this.flags[i])) { changed = true; }
			if (this.checkboxes[i].checked) {
				settings.flags |= this.flags[i];
			} else {
				settings.flags &= ~this.flags[i];
			}
		}
		return changed;
	}
	loadFromSettings(settings) {
		let changed = false;
		for (let i = 0; i < this.checkboxes.length; i++) {
			if (this.checkboxes[i].checked != !!(settings.flags & this.flags[i])) { changed = true; }
			this.checkboxes[i].checked = !!(settings.flags & this.flags[i]);
		}
		return changed;
	}
}

class RadioButtonsControl {
	constructor(parentElement, flagBit, flagMask, name, labels) {
		this.flagBit = flagBit;
		this.flagMask = flagMask;
		this.onchange = null;

		this.rootElement = addElement(parentElement, 'div', 'control-input-bag');
		this.buttons = [];
		labels.forEach(labelRow => {
			let row = addElement(this.rootElement, 'div', 'control-input-bag-row');
			labelRow.forEach(label => {
				let labelNode = addLabel(row, 'control-input-bag-contents', label);
				let button = addInput(labelNode, '', 'radio');
				button.name = name;
				button.value = this.buttons.length;
				button.addEventListener('change', () => {
					if (this.onchange) { this.onchange(button.value); }
				});
				this.buttons.push(button);
			});
		});
	}
	storeState() {
		for (let i = 0; i < this.buttons.length; i++) {
			if (this.buttons[i].checked) { return i; }
		}
		return 0;
	}
	loadState(state) {
		this.buttons[state].checked = true;
	}
	storeToSettings(settings) {
		let setting = (settings.flags >> this.flagBit) & this.flagMask;
		let changed = !this.buttons[setting].checked;
		for (let i = 0; i < this.buttons.length; i++) {
			if (this.buttons[i].checked) {
				settings.flags &= ~(this.flagMask << this.flagBit);
				settings.flags |= i << this.flagBit;
				break;
			}
		}
		return changed;
	}
	loadFromSettings(settings) {
		let setting = (settings.flags >> this.flagBit) & this.flagMask;
		let changed = !this.buttons[setting].checked;
		this.buttons[setting].checked = true;
		return changed;
	}
}

class RangeWithTextFieldControl {
	constructor(parentElement, setting, min, max, step, rangeToTextValueFn) {
		this.setting = setting;
		this.value = 0.0;
		this.onchange = null;

		this.rootElement = addElement(parentElement, 'div', 'control-input-slider');
		let label0 = addLabel(this.rootElement, 'control-input-slider', '');
		this.range = addInput(label0, 'control-range', 'range');
		this.range.min = min;
		this.range.max = max;
		this.range.step = step || 'any';
		this.range.value = this.value;
		let label1 = addLabel(this.rootElement, '', '');
		this.text = addInput(label1, 'control-range-text', 'text');
		this.text.value = rangeToTextValueFn(this.value);

		this.formatValueForText = (value) => {
			let str = ''+parseFloat(value);
			if (str.match('e')) { str = parseFloat(value).toFixed(20); } // No scientific notation allowed!
			let strBuilder = [];
			let firstSigIdx = 100000;
			let decimalIdx = 100000;
			for (let i = 0; i < str.length; i++) {
				if (i < firstSigIdx && str[i] >= '1' && str[i] <= '9') { firstSigIdx = i; }
				if (i < decimalIdx && str[i] == '.') { decimalIdx = i; ++firstSigIdx; }
				if (i >= decimalIdx && i > firstSigIdx + 2) { break; }
				strBuilder.push(i < firstSigIdx + 3 ? str[i] : '0');
			}
			while (strBuilder.length > decimalIdx && strBuilder[strBuilder.length - 1] == '0') { strBuilder.pop(); }
			if (strBuilder[strBuilder.length - 1] == '.') { strBuilder.pop(); }
			return strBuilder.join('');
		};
		this.searchValueForRange = (value) => {
			let low = min;
			let high = max;
			let guess = 0.5 * (low + high);
			for (let itr = 0; itr < 32; itr++) {
				if (value < rangeToTextValueFn(guess)) {
					high = guess;
				} else {
					low = guess;
				}
				guess = 0.5 * (low + high);
			}
			return guess;
		};
		this.range.addEventListener('input', () => {
			let value = rangeToTextValueFn(this.range.value);
			this.text.value = this.formatValueForText(value);
			this.value = parseFloat(this.text.value);
			if (this.onchange) { this.onchange(this.value); }
		});
		this.text.addEventListener('change', () => {
			let isNumeric = !isNaN(this.text.value) && !isNaN(parseFloat(this.text.value));
			if (!isNumeric) { return; }
			this.value = parseFloat(this.text.value);
			this.range.value = this.searchValueForRange(this.value);
			if (this.onchange) { this.onchange(this.value); }
		});
	}
	storeState() {
		return { value: this.value, range: this.range.value, text: this.text.value };
	}
	loadState(state) {
		this.value = state.value;
		this.range.value = state.range;
		this.text.value = state.text;
	}
	storeToSettings(settings) {
		let changed = settings[this.setting] != this.value;
		settings[this.setting] = this.value;
		return changed;
	}
	loadFromSettings(settings) {
		let changed = this.value != settings[this.setting];
		this.value = settings[this.setting];
		this.text.value = this.formatValueForText(this.value);
		this.range.value = this.searchValueForRange(this.value);
		return changed;
	}
}

function addControl(parentElement, label, createControlFn) {
	let row = addElement(parentElement, 'div', 'control-row');
	let labelNode = addElement(row, 'div', 'control-label');
	labelNode.innerHTML = label;

	let inputArea = addElement(row, 'div', 'control-input-area');

	let ctrl0 = createControlFn(inputArea, 0);
	ctrl0.loadFromSettings(settings);

	let link = addElement(inputArea, 'div', 'control-link');
	let linkLabel = addLabel(link, 'link-toggle', '');
	let linkCheckbox = addInput(linkLabel, 'link-toggle', 'checkbox');
	let linkText = addElement(linkLabel, 'span', 'link-text');
	links.push({linkCheckbox, label}); // Hackily push the link into the global link list

	let ctrl1 = createControlFn(inputArea, 1);
	ctrl1.loadFromSettings(settings1);

	let undoState = null;

	ctrl0.onchange = (value) => {
		let changed = ctrl0.storeToSettings(settings);
		if (changed) { updateSettings(0); console.log(0, label, ctrl0.storeState()); }
		if (linkCheckbox.checked) {
			ctrl1.loadState(ctrl0.storeState());
			undoState = ctrl1.storeState();
			let changed1 = ctrl1.storeToSettings(settings1);
			if (changed1) { updateSettings(1); console.log(1, label, ctrl1.storeState()); }
		}
	};
	ctrl1.onchange = (value) => {
		undoState = ctrl1.storeState();
		let changed = ctrl1.storeToSettings(settings1);
		if (changed) { updateSettings(1); console.log(1, label, ctrl1.storeState()); }
		// Disable the link
		linkCheckbox.checked = false;
		ctrl1.rootElement.classList.remove('control-input-linked');
	};

	// Set up the link checkbox functionality
	linkCheckbox.addEventListener('change', () => {
		if (linkCheckbox.checked) {
			ctrl1.rootElement.classList.add('control-input-linked');
			ctrl1.loadState(ctrl0.storeState());
		} else {
			ctrl1.rootElement.classList.remove('control-input-linked');
			if (undoState) { ctrl1.loadState(undoState); }
		}
		let changed1 = ctrl1.storeToSettings(settings1);
		if (changed1) { updateSettings(1); console.log(1, label, 'link', ctrl1.storeState()); }
	});
}

function addCheckbox(parentElement, label, flags, checkboxLabels) {
	addControl(parentElement, label, (parent, side) => {
		return new CheckboxControl(parent, flags, checkboxLabels);
	});
}

function addRadioButtons(parentElement, label, flagBit, flagMask, buttonLabels) {
	addControl(parentElement, label, (parent, side) => {
		return new RadioButtonsControl(parent, flagBit, flagMask, label+side, buttonLabels);
	});
}

function addRangeWithTextField(parentElement, label, setting, min, max, step, rangeToTextValueFn) {
	addControl(parentElement, label, (parent, side) => {
		return new RangeWithTextFieldControl(parent, setting, min, max, step, rangeToTextValueFn);
	});
}

document.addEventListener('DOMContentLoaded', function(event) {
	// Gather all the slide dom elements
	// slidesChildren = document.getElementById('slides').childNodes;
	// for (let i = 0; i < slidesChildren.length; i++) {
	// 	if (slidesChildren[i].className == 'slide') {
	// 		slides.push(slidesChildren[i]);
	// 	}
	// }

	// // Add table of contents entries for all the slides
	// let toc = document.getElementById('toc');
	// function addTocEntry(label, onclick) {
	// 	let span = document.createElement('span');
	// 	let text = document.createTextNode(label);
	// 	span.appendChild(text);
	// 	toc.appendChild(span);
	// 	span.className = 'toc-entry';
	// 	span.addEventListener('click', onclick);
	// }
	// addTocEntry('<', () => gotoSlide(currentSlide - 1));
	// for (let i = 0; i < slides.length; i++) {
	// 	addTocEntry(''+(i + 1), ((idx) => () => gotoSlide(idx))(i));
	// }
	// addTocEntry('>', () => gotoSlide(currentSlide + 1));

	// Support navigating through slides with the keyboard
	document.addEventListener('keydown', function(event) {
		// We never want to look at input if the user is typing in a text box
		if (document.activeElement.tagName.toLowerCase() == 'input' &&
			document.activeElement.type.toLowerCase() == 'text') { return; }

		if (event.key.toUpperCase() == 'R') {
			resetBlocks();
		}
		if (event.key.toUpperCase() == 'L') {
			toggleAllLinks();
		}

		// Furthermore, certain inputs we don't want to look at if the user has anything selected
		if (document.activeElement != document.body) { return; }

		// if (event.key == 'ArrowLeft' || event.key == 'ArrowRight') {
		// 	if (event.key == 'ArrowLeft') {
		// 		gotoSlide(currentSlide - 1);
		// 	} else if (event.key == 'ArrowRight') {
		// 		gotoSlide(currentSlide + 1);
		// 	}
		// 	event.preventDefault();
		// }
	}, true);

	//
	// gotoSlide(0);
});
