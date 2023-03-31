import json
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def json_to_pdf(data, output_path):
    fields = data["features"][0]["properties"]
    c = canvas.Canvas(output_path, pagesize=letter)
    text = c.beginText()
    text.setFont("Helvetica", 12)
    text.setTextOrigin(50, 750)
    text.textLine("Shapefile Field Information:")
    for field_name, field_info in fields.items():
        field_type = field_info["type"]
        text.textLine(f"{field_name}: {field_type}")
    c.drawText(text)
    c.save()

try:
    with open('metadata.json', 'r') as f:
        data_str = f.read()
        data = json.loads(json.loads(data_str))
        # data = eval(data)
        print(data)
except FileNotFoundError:
    print("metadata.json file not found")
    exit()
except json.JSONDecodeError as e:
    print("Error parsing JSON:", e)
    exit()

if not isinstance(data, dict):
    print("Loaded data is not a dictionary")
    exit()

json_to_pdf(data, 'metadata_pdf.pdf')
