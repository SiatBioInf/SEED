import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
from openai import OpenAI

client = OpenAI(
    api_key=os.environ["SEED_LLM_API_KEY"],
    base_url="https://uni-api.cstcloud.cn/v1"
)

# Simple test
resp = client.chat.completions.create(
    model="S1-Base-Ultra",
    messages=[{"role": "user", "content": "Say 'API works' in Chinese"}],
    temperature=0,
    max_tokens=50,
    timeout=30
)
print(resp.choices[0].message.content)
